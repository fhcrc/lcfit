#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/Io/IoSequenceFactory.h>
#include <Bpp/Seq/SiteTools.h>

#include "lcfit.h"
#include "gsl.h"

using namespace std;

const size_t MAX_EVAL_COUNT = 100;
typedef bpp::TreeTemplate<bpp::Node> Tree;

template<typename T>
void print_vector(std::vector<T> v, const char delim = '\t', ostream& out = std::cout)
{
    for(auto & i : v) {
        out << i << delim;
    }
    out << endl;
}

// A simple point
struct Point
{
    Point(const double x, const double y) :
        x(x),
        y(y) {};
    double x, y;
};

// Support hashing points
namespace std {
    template <> struct hash<Point> {
        size_t operator()(const Point & p) const
        {
            return std::hash<double>()(p.x) ^ std::hash<double>()(p.y);
        };
    };
}

/// Evaluation of the fit for a node at a given branch length
struct Evaluation {
    Evaluation(int node_id, double branch_length, double ll, double pred_ll) :
        node_id(node_id),
        branch_length(branch_length),
        ll(ll),
        pred_ll(pred_ll) {};
    int node_id;
    double branch_length, ll, pred_ll;
};

void to_csv(ostream& out, const Evaluation& e)
{
    out << e.node_id << ","
        << e.branch_length << ","
        << e.ll << ","
        << e.pred_ll << endl;
}

enum class Monotonicity
{
    UNKNOWN = 0,
    MONO_INC,
    MONO_DEC,
    NON_MONOTONIC
};

/// Find the monotonicity of a set of (sorted) points
inline Monotonicity monotonicity(const vector<Point>& points)
{
    assert(points.size() > 0);
    auto i = begin(points);
    auto e = end(points);
    bool maybe_inc = true, maybe_dec = true;

    Point last = *i++;
    for(; i != e; i++) {
        Point current = *i;
        if(current.y > last.y) maybe_dec = false;
        else if(current.y < last.y) maybe_inc = false;
        last = current;
    }
    assert(!(maybe_inc && maybe_dec));

    if(!maybe_inc && !maybe_dec) return Monotonicity::NON_MONOTONIC;
    else if(maybe_inc) return Monotonicity::MONO_INC;
    else if(maybe_dec) return Monotonicity::MONO_DEC;
    assert(false);
}

template<typename ForwardIterator>
inline size_t max_index(ForwardIterator first, const ForwardIterator last)
{
    assert(first != last);
    ForwardIterator max_iter = first++;
    size_t max_idx = 0;
    for(size_t idx = 0; first != last; ++first, ++idx) {
        if(*first > *max_iter) {
            max_iter = first;
            max_idx = idx;
        }
    }

    return max_idx;
}

template<typename ForwardIterator>
inline size_t min_index(ForwardIterator first, const ForwardIterator last)
{
    assert(first != last);
    ForwardIterator min_iter = first++;
    size_t min_idx = 0;
    for(size_t idx = 0; first != last; ++first, ++idx) {
        if(*first < *min_iter) {
            min_iter = first;
            min_idx = idx;
        }
    }

    return min_idx;
}

// Actual comparison code

/// Calculates the log-likelihood of a tree given data and a model
class TreeLikelihoodCalculator
{
public:
    TreeLikelihoodCalculator(bpp::SiteContainer* sites, bpp::SubstitutionModel* model, bpp::DiscreteDistribution* rate_dist) :
        sites(sites),
        model(model),
        rate_dist(rate_dist) {};

    /// Calculate log-likelihood
    double calculate_log_likelihood(const Tree& tree) const
    {
        // First, sort out the appropriate likelihood calculator
        // Approximates bppml behavior
        unique_ptr<bpp::DiscreteRatesAcrossSitesTreeLikelihood> calc;
        if(dynamic_cast<bpp::MixedSubstitutionModel*>(model) == 0) {
            calc.reset(
                    new bpp::RHomogeneousTreeLikelihood(tree, *sites, model, rate_dist, false, false, false));
        } else {
            // Mixed model
            calc.reset(
                    new bpp::RHomogeneousMixedTreeLikelihood(tree, *sites, model, rate_dist, false, false, false));
        }
        calc->initialize();
        //like->computeTreeLikelihood();
        return calc->getLogLikelihood();
    }
private:
    bpp::SiteContainer* sites;
    bpp::SubstitutionModel* model;
    bpp::DiscreteDistribution* rate_dist;
};

struct LCFitResult
{
    vector<double> coef;
    vector<Point> evaluated_points;
    size_t ll_eval_count;
};

/// Runs lcfit, generating parameters for the CFN model
class LCFitter
{
public:
    LCFitter(const vector<double> start, const vector<double> sample_points, TreeLikelihoodCalculator* calc, ostream* csv_fit_out) :
        start(start),
        sample_points(sample_points),
        calc(calc),
        csv_fit_out(csv_fit_out) {};

    /// Run lcfit, returning coefficients of the model
    ///
    /// Given a set of starting points, attempts to add (branch_length, ll) points until the function is non-monotonic,
    /// then fits the CFN model using these sampled points.
    LCFitResult fit_model(const Tree& tree, const int node_id, FILE* log_fp=nullptr)
    {
        vector<double> x = start; // Initial conditions for [c,m,r,b]
        const vector<Point> points = this->select_points(tree, node_id);

        // Log fit
        if(csv_fit_out != nullptr) {
            for(const Point& p : points) {
                *csv_fit_out << node_id << "," << p.x << "," << p.y << endl;
            }
        }

        // Scale initial conditions to intersect with maximum likelihood point
        scale_coefs(x, points);

        // Extract x, y
        vector<double> t, l;
        t.reserve(points.size());
        std::transform(begin(points), end(points), std::back_inserter(t), [](const Point& p) ->double { return p.x; });
        l.reserve(points.size());
        std::transform(begin(points), end(points), std::back_inserter(l), [](const Point& p) ->double { return p.y; });

        const int status = fit_ll_log(t.size(), t.data(), l.data(), x.data(), log_fp);
        if(status) throw runtime_error("fit_ll returned: " + std::to_string(status));
        return {x, points, points.size()};
    }

    pair<double, size_t> estimate_ml_branch_length(Tree tree, const int node_id, const double tol)
    {
        vector<double> x = start; // Initial conditions for [c,m,r,b]
        const vector<Point> points = this->select_points(tree, node_id);
        scale_coefs(x, points);
        size_t eval_count = points.size();


        // Extract x, y
        vector<double> t, l;
        t.reserve(points.size());
        std::transform(begin(points), end(points), std::back_inserter(t), [](const Point& p) ->double { return p.x; });
        l.reserve(points.size());
        std::transform(begin(points), end(points), std::back_inserter(l), [](const Point& p) ->double { return p.y; });
        double last = ml_t(x[0], x[1], x[2], x[3]);
        while(eval_count <= MAX_EVAL_COUNT) {
            const int status = fit_ll(t.size(), t.data(), l.data(), x.data());
            if(status) throw runtime_error("fit_ll returned: " + std::to_string(status));
            double ml_bl = ml_t(x[0], x[1], x[2], x[3]);

            // Return minimum branch length
            if(ml_bl < 0) return {1e-7, eval_count};

            // Close enough?
            if(std::abs(ml_bl - last) < tol) return {ml_bl, eval_count};

            //cerr << node_id << ": " << ml_bl << endl;

            // Add current ML distance to fit
            t.push_back(ml_bl);
            tree.setDistanceToFather(node_id, ml_bl);
            l.push_back(calc->calculate_log_likelihood(tree));
            eval_count++;
            // Subset to top 4 LL values
            this->keep_top(t, l, 4); // TODO: Fix magic number 4?
            eval_count++;
            last = ml_bl;
        }
        return {last, eval_count};
    };

private:
    const vector<double> start;
    const vector<double> sample_points;
    TreeLikelihoodCalculator* calc;
    ostream* csv_fit_out;

    // Keep top `n` y values
    void keep_top(vector<double>& x, vector<double>& y, const size_t n)
    {
        if(n >= y.size()) return;
        auto xbeg = begin(x), xend = end(x);
        auto ybeg = begin(y), yend = end(y);
        vector<double> sorted_y = y;
        std::sort(begin(sorted_y), end(sorted_y));
        const double y_cutoff = sorted_y[sorted_y.size() - n];
        for(; xbeg != xend && ybeg != yend; ++xbeg, ++ybeg) {
            if(*ybeg < y_cutoff) {
                y.erase(ybeg);
                x.erase(xbeg);
            }
        }
    };

    // Scale initial conditions to intersect with maximum likelihood point
    void scale_coefs(vector<double>& x, const vector<Point>& points)
    {
        auto p = std::max_element(begin(points), end(points),
                [](const Point& p1, const Point& p2) -> bool { return p1.y > p2.y; });
        double scale_factor = cm_scale_factor(p->x, p->y, x[0], x[1], x[2], x[3]);
        x[0] *= scale_factor;
        x[1] *= scale_factor;
    };

    /// Choose the input (branch_length, ll) samples for running lcfit
    vector<Point> select_points(Tree tree, const int node_id, const size_t max_points=8) // TODO: Fix magic number 8
    {
        vector<Point> points;

        // Try starting points
        for(const double & d : sample_points) {
            assert(d >= 0);
            tree.setDistanceToFather(node_id, d);
            points.emplace_back(d, calc->calculate_log_likelihood(tree));
        }

        // Add additional samples until the evaluated branch lengths enclose a maximum.
        size_t offset; // Position - points kept sorted by x-value
        double d;      // Branch length

        Monotonicity c = monotonicity(points);
        do {
            switch(c) {
            case Monotonicity::NON_MONOTONIC:
                d = (points[1].x + points[2].x) / 2.0; // Add a point between the first and second try
                offset = 2;
                break;
            case Monotonicity::MONO_INC:
                d = points.back().x * 2.0; // Double largest value
                offset = points.size();
                break;
            case Monotonicity::MONO_DEC:
                d = points[0].x / 10.0; // Add new smallest value - order of magnitude lower
                offset = 0;
                break;
            default:
                assert(false);
            }

            assert(d >= 0);
            tree.setDistanceToFather(node_id, d);
            points.emplace(points.begin() + offset, d, calc->calculate_log_likelihood(tree));

            c = monotonicity(points);

            assert(is_sorted(points.begin(), points.end(),
                        [](const Point& p1, const Point& p2) -> bool { return p1.x < p1.y; }));
        } while(points.size() <= max_points && c != Monotonicity::NON_MONOTONIC);

        return points;
    };
};

// Evaluate log-likelihood obtained by model fit versus actual log-likelihood at a variety of points
vector<Evaluation> evaluate_fit(Tree tree, TreeLikelihoodCalculator* calc, const int node_id, const vector<double>& x, const double max_t=1.0, const double delta=0.01) {
    vector<Evaluation> evaluations;
    const double c = x[0];
    const double m = x[1];
    const double r = x[2];
    const double b = x[3];
    const double t = tree.getDistanceToFather(node_id);
    for(double t = delta; t <= max_t; t += delta) {
        tree.setDistanceToFather(node_id, t);
        double actual_ll = calc->calculate_log_likelihood(tree);
        double fit_ll = ll(t, c, m, r, b);
        evaluations.emplace_back(node_id, t, actual_ll, fit_ll);
    }
    tree.setDistanceToFather(node_id, t);
    return evaluations;
}

pair<double, size_t> estimate_ml_branch_length_brent(const TreeLikelihoodCalculator& calc,
        Tree tree,
        const size_t node_id,
        double left,
        double right,
        double raw_start,
        double tolerance) {
    size_t n_eval=0;

    std::function<double(double)> f = [&node_id, &calc, &tree, &n_eval](double d)->double {
        n_eval++;
        assert(d >= 0);
        tree.setDistanceToFather(node_id, d);

        // Minimize -ll
        return -calc.calculate_log_likelihood(tree);
    };

    double lefty = f(left), righty = f(right);
    double miny = std::min(lefty, righty);
    double smaller = lefty < righty ? left : right;

    assert(left >= 0);
    assert(right >= 0);

    double prev_start = raw_start;
    for(size_t iteration = 0; ; iteration++) {
        assert(prev_start >= 0);
        if(std::abs(prev_start - smaller) < tolerance)
            return {smaller, n_eval}; // Abutting `left`

        double prev_val = f(prev_start);
        if(prev_val < miny) {
            break; /* Appropriate starting point */
        } else if(iteration > MAX_EVAL_COUNT) {
            throw std::runtime_error("Minimization_brent exceeded max iters");
        } else {
            prev_start = (prev_start + smaller) / 2;
        }
    }

    double ml_bl = gsl::minimize(f, prev_start, left, right);

    // gsl::minimize calls f(left), f(right), and f(prev_start)
    // If we were smarter, could cache those.
    // To make it fair, correct for double counting.
    return {ml_bl, n_eval - 3};
}

int run_main(int argc, char** argv)
{
    bpp::BppApplication lcfit_compare(argc, argv, "lcfit-compare");

    map<string, string> params = lcfit_compare.getParams();

    /********************/
    /* ARGUMENT PARSING */
    /********************/

    // Alphabet
    unique_ptr<bpp::Alphabet> alphabet(bpp::SequenceApplicationTools::getAlphabet(params, "", false));
    // Sites
    unique_ptr<bpp::VectorSiteContainer> all_sites(bpp::SequenceApplicationTools::getSiteContainer(alphabet.get(), params));
    unique_ptr<bpp::VectorSiteContainer> sites(bpp::SequenceApplicationTools::getSitesToAnalyse(*all_sites, params, "", true, false));
    all_sites.reset();
    // Tree
    unique_ptr<bpp::Tree> in_tree(bpp::PhylogeneticsApplicationTools::getTree(params));
    bpp::TreeTemplate<bpp::Node> tree(*in_tree);
    // Model
    unique_ptr<bpp::SubstitutionModel> model(bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet.get(), sites.get(), params));
    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);
    // Rate dist
    unique_ptr<bpp::DiscreteDistribution> rate_dist(bpp::PhylogeneticsApplicationTools::getRateDistribution(params));


    double fit_tolerance = bpp::ApplicationTools::getDoubleParameter("lcfit.ml_tolerance", params, 1e-4);
    // lcfit-specific
    // Output
    string csv_like_path = bpp::ApplicationTools::getAFilePath("lcfit.output.likelihoods_file", params, true, false);
    ofstream csv_like_out(csv_like_path);
    string csv_ml_path = bpp::ApplicationTools::getAFilePath("lcfit.output.maxima_file", params, true, false);
    ofstream csv_ml_out(csv_ml_path);
    string csv_fit_path = bpp::ApplicationTools::getAFilePath("lcfit.output.fit_file", params, true, false);
    ofstream csv_fit_out(csv_fit_path);
    csv_fit_out << "node_id,branch_length,ll" << endl;;

    FILE* fit_log_fp = nullptr;
    if(bpp::ApplicationTools::parameterExists("lcfit.output.fit_log", params)) {
        fit_log_fp = fopen(bpp::ApplicationTools::getAFilePath("lcfit.output.fit_log", params, true, false).c_str(), "w");
        assert(fit_log_fp != nullptr);
        fprintf(fit_log_fp, "iter,c,m,r,b,status\n");
    }

    const vector<double> sample_points = bpp::ApplicationTools::getVectorParameter<double>(
            "lcfit.sample.branch.lengths",
            params, ',', "0.1,0.15,0.5");
    const vector<double> start = bpp::ApplicationTools::getVectorParameter<double>("lcfit.starting.values",
                                 params, ',', "1500,1000,2.0,0.5");

    // Calculators
    TreeLikelihoodCalculator likelihood_calc(sites.get(), model.get(), rate_dist.get());
    LCFitter fit(start, sample_points, &likelihood_calc, &csv_fit_out);

    /*
     * Run evaluations on each node, write output
     */
    csv_like_out << "node_id,branch_length,bpp_ll,fit_ll" << endl;
    //csv_ml_out << "node_id,t,t_hat,c,m,r,b" << endl;
    csv_ml_out << "node_id,lcfit_t,lcfit_n,brent_t,brent_n,lcfit_fit_t,lcfit_fit_n" << endl;
    for(const int & node_id : tree.getNodesId()) {
        cerr << "Node " << node_id << "\r";
        if(!tree.hasDistanceToFather(node_id)) continue;
        LCFitResult r = fit.fit_model(tree, node_id, fit_log_fp);

        // Find maximum of all x-values and 1.0
        double max_x = accumulate(begin(r.evaluated_points),
                                  end(r.evaluated_points),
                                  1.0,
                                  [](const double m, const Point& p) { return std::max(m, p.x); });

        const vector<Evaluation> evals = evaluate_fit(tree, &likelihood_calc, node_id, r.coef, max_x,
                                                      max_x/100.0);
        // Write to CSV
        std::for_each(begin(evals), end(evals),
            [&csv_like_out](const Evaluation & e) { to_csv(csv_like_out, e); });

        auto lcfit_ml_result = fit.estimate_ml_branch_length(tree, node_id, fit_tolerance);
        auto brent_ml_result = estimate_ml_branch_length_brent(likelihood_calc,
                tree,
                node_id,
                1e-6,
                3.,
                0.5, fit_tolerance);
        csv_ml_out << node_id << ","
            << lcfit_ml_result.first << ","
            << lcfit_ml_result.second << ","
            << brent_ml_result.first << ","
            << brent_ml_result.second << ","
            << ml_t(r.coef[0], r.coef[1], r.coef[2], r.coef[3]) << "," <<
            r.evaluated_points.size() << endl;
    }

    return 0;
}

int main(int argc, char** argv)
{
    try {
        run_main(argc, argv);
    } catch(exception& e) {
        cerr << "Error:" << e.what() << endl;
        return 1;
    }
}
