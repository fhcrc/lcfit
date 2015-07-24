#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
#include <utility>

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Io/IoSequenceFactory.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/SiteTools.h>

#include "gsl.h"
#include "lcfit.h"
#include "lcfit_select.h"
#include "lcfit_cpp.h"

using namespace std;
using namespace lcfit;

typedef bpp::TreeTemplate<bpp::Node> Tree;

template<typename T>
void print_vector(std::vector<T> v, const char delim = '\t', ostream& out = std::cout)
{
    for(auto & i : v) {
        out << i << delim;
    }
    out << endl;
}


/// LogLikelihoodComparison of the fit for a node at a given branch length
struct LogLikelihoodComparison {
    LogLikelihoodComparison(int node_id, double branch_length, double ll, double pred_ll) :
        node_id(node_id),
        branch_length(branch_length),
        ll(ll),
        pred_ll(pred_ll) {};
    int node_id;
    double branch_length, ll, pred_ll;
};

/// Model fit
struct ModelFit {
    ModelFit(int node_id, double t, double t_hat, double c, double m, double r, double b) :
        node_id(node_id),
        t(t),
        t_hat(t_hat),
        c(c),
        m(m),
        r(r),
        b(b) {};
    int node_id;
    double t,     // Branch length
           t_hat, // Fit branch length
           c,
           m,
           r,
           b;
};

void to_csv(ostream& out, const LogLikelihoodComparison& e)
{
    out << e.node_id << ","
        << e.branch_length << ","
        << e.ll << ","
        << e.pred_ll << endl;
}

void to_csv(ostream& out, const ModelFit& f)
{
    out << f.node_id << ","
        << f.t << ","
        << f.t_hat << ","
        << f.c << ","
        << f.m << ","
        << f.r << ","
        << f.b << endl;
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
    TreeLikelihoodCalculator(const bpp::TreeLikelihood& tl) :
        like(tl.clone())
    {
        like->initialize();
    };

    /// Calculate log-likelihood
    double calculate_log_likelihood()
    {
        //like->computeTreeLikelihood();
        return like->getLogLikelihood();
    }

    void set_branch_length(const size_t node, double length)
    {
        length = std::max(length, 1e-6);
        like->setParameterValue("BrLen" + std::to_string(node), length);
    }

    void get_branch_length(const size_t node) const
    {
        like->getParameterValue("BrLen" + std::to_string(node));
    }
//private:
    std::unique_ptr<bpp::TreeLikelihood> like;
};

/// Calculate the log likelihood of a node given a branch length
struct NodeLikelihoodCalculator
{
    NodeLikelihoodCalculator(TreeLikelihoodCalculator *c,
                             const int node_id) :
        calc(c),
        node_id(node_id),
        n_evals(0)
    {};

    TreeLikelihoodCalculator* calc;
    const int node_id;
    size_t n_evals;

    double operator()(double d)
    {
        n_evals++;
        calc->set_branch_length(node_id, d);
        return calc->calculate_log_likelihood();
    };

};

double lcfit_ll(double t, void *data)
{
    NodeLikelihoodCalculator* calc = static_cast<NodeLikelihoodCalculator*>(data);
    return (*calc)(t);
}

/// Runs lcfit, generating parameters for the BSM
class LCFitter
{
public:
    LCFitter(const vector<double>& start, const vector<double>& sample_points, TreeLikelihoodCalculator* calc, ostream* csv_fit_out) :
        start(start),
        sample_points(sample_points),
        calc(calc),
        csv_fit_out(csv_fit_out) {};


    /// Run lcfit, returning coefficients of the model
    ///
    /// Given a set of starting points, attempts to add (branch_length, ll) points until the function is non-monotonic,
    /// then fits the BS model using these sampled points.
    bsm_t fit_model(const int node_id, const int max_iter = 250)
    {
        bsm_t m = {start[0], start[1], start[2], start[3]};
        NodeLikelihoodCalculator ll(calc,node_id);
        lcfit::LCFitResult fit_result = lcfit::fit_bsm_log_likelihood(ll, m, sample_points, 8, max_iter);

        // Log fit
        if(csv_fit_out != nullptr) {
            for(const Point& p : fit_result.evaluated_points) {
                *csv_fit_out << node_id << "," << std::setprecision(20) << p.x << "," << p.y << endl;
            }
        }

        // Scale initial conditions to intersect with maximum likelihood point
        return fit_result.model_fit;
    }

    pair<double, size_t> estimate_ml_branch_length(const int node_id, const double tol=1e-5)
    {
        bsm_t m = {start[0], start[1], start[2], start[3]};
        NodeLikelihoodCalculator ll(calc, node_id);
        log_like_function_t fn { lcfit_ll, &ll };
        bool success = 1;

        /* GSL requires at least four points for minimization. */
        std::vector<Point> evaluated_points = lcfit::select_points(ll, sample_points, 4);
        std::vector<double> pts;
        for (const auto& point : evaluated_points) {
            pts.push_back(point.x);
        }

        const double result = estimate_ml_t(&fn,
                                            pts.data(),
                                            pts.size(),
                                            tol,
                                            &m,
                                            &success);

        return pair<double, size_t>(result, ll.n_evals);
    }

    pair<double, size_t> estimate_ml_branch_length_brent(const int node_id,
                                                         const double tolerance=1e-5,
                                                         double left=1e-6,
                                                         double right=2,
                                                         double raw_start=0.1)
    {
        size_t n_eval = 0;

        NodeLikelihoodCalculator ll(calc, node_id);
        std::function<double(double)> f = [&ll, &n_eval](double d)->double {
            n_eval++;
            return -ll(d);
        };

        double lefty = f(left), righty = f(right);
        double miny = std::min(lefty, righty);
        double smaller = lefty < righty ? left : right;

        double prev_start = raw_start;
        for(size_t iteration = 0; ; iteration++) {

            if(std::abs(prev_start - smaller) < tolerance)
                return {smaller, n_eval}; // Abutting `left`

            double prev_val = f(prev_start);
            if(prev_val < miny) {
                break; /* Appropriate starting point */
            } else if(iteration > 100) {
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
private:
    const vector<double> start;
    const vector<double> sample_points;
    TreeLikelihoodCalculator* calc;
    ostream* csv_fit_out;
};

// Evaluate log-likelihood obtained by model fit versus actual log-likelihood at a variety of points
vector<LogLikelihoodComparison> compare_log_likelihoods(Tree tree, TreeLikelihoodCalculator* calc, const int node_id, const bsm_t& m)
{
    vector<LogLikelihoodComparison> evaluations;
    const double ml_t = tree.getDistanceToFather(node_id);
    const double ml_ll = calc->calculate_log_likelihood();
    const double threshold = ml_ll - 100.0;

    // Note this function's side effect of setting node_id's branch
    // length to t!
    //
    // Also, would the threshold be better as a ratio of ml_ll instead
    // of a constant offset?
    auto bounds_fn = [&calc, node_id, threshold](double t) {
        calc->set_branch_length(node_id, t);
        return calc->calculate_log_likelihood() - threshold;
    };

    const int MAX_ITER = 100;
    const double TOLERANCE = 1e-3;

    // Limits are from Bio++: minimum branch length which may be considered is 1e-6
    double lower = 1e-6;

    // refine lower bound if possible
    if (ml_t > lower && bounds_fn(lower) < 0.0) {
        lower = gsl::find_root(bounds_fn, lower, ml_t, MAX_ITER, TOLERANCE);
    }

    double upper_min = ml_t;
    double upper_max = ml_t * 2.0;

    while (bounds_fn(upper_max) > 0.0) {
        upper_min = upper_max;
        upper_max *= 2.0;
    }

    double upper = gsl::find_root(bounds_fn, upper_min, upper_max, MAX_ITER, TOLERANCE);

    const size_t N_SAMPLES = 501;
    const double delta = (upper - lower) / (N_SAMPLES - 1);

    for (size_t i = 0; i < N_SAMPLES; ++i) {
        const double t = lower + (i * delta);
        calc->set_branch_length(node_id, t);
        double actual_ll = calc->calculate_log_likelihood();
        double fit_ll = lcfit_bsm_log_like(t, &m);
        evaluations.emplace_back(node_id, t, actual_ll, fit_ll);
    }
    calc->set_branch_length(node_id, ml_t);
    assert(tree.getDistanceToFather(node_id) == ml_t);
    return evaluations;
}

// Compare ML branch length estimate from lcfit to original branch length (presumed to be ML value)
ModelFit compare_ml_values(const Tree& tree, const int node_id, const bsm_t m)
{
    const double t_hat = lcfit_bsm_ml_t(&m);
    const double t = tree.getDistanceToFather(node_id);
    return ModelFit(node_id, t, t_hat, m.c, m.m, m.r, m.b);
}

int run_main(int argc, char** argv)
{
    bpp::BppApplication lcfit_compare(argc, argv, "lcfit-compare");

    map<string, string> params = lcfit_compare.getParams();

    /********************/
    /* ARGUMENT PARSING */
    /********************/

    // Alphabet
    const unique_ptr<bpp::Alphabet> alphabet(bpp::SequenceApplicationTools::getAlphabet(params, "", false));
    // Genetic code
    unique_ptr<bpp::GeneticCode> gcode;
    if(bpp::CodonAlphabet* codon_alph = dynamic_cast<bpp::CodonAlphabet*>(alphabet.get())) {
        const string code_desc = bpp::ApplicationTools::getStringParameter("genetic_code", params, "Standard", "", true, true);
        bpp::ApplicationTools::displayResult("Genetic Code", code_desc);
        gcode.reset(bpp::SequenceApplicationTools::getGeneticCode(codon_alph->getNucleicAlphabet(), code_desc));
    }
    // Sites
    unique_ptr<bpp::VectorSiteContainer> all_sites(bpp::SequenceApplicationTools::getSiteContainer(alphabet.get(), params));
    unique_ptr<bpp::VectorSiteContainer> sites(bpp::SequenceApplicationTools::getSitesToAnalyse(*all_sites, params, "", true, false));
    all_sites.reset();
    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    // Tree
    unique_ptr<bpp::Tree> in_tree(bpp::PhylogeneticsApplicationTools::getTree(params));
    bpp::TreeTemplate<bpp::Node> tree(*in_tree);
    // Rate dist
    unique_ptr<bpp::DiscreteDistribution> rate_dist(bpp::PhylogeneticsApplicationTools::getRateDistribution(params));
    // Model
    // See bppSeqGen.cpp L253
    const std::string nh_opt = bpp::ApplicationTools::getStringParameter("nonhomogeneous", params, "no", "", true, false);
    unique_ptr<bpp::SubstitutionModel> model;
    unique_ptr<bpp::SubstitutionModelSet> model_set;
    unique_ptr<bpp::TreeLikelihood> bpp_tree_like;
    if(nh_opt == "no") {
        model.reset(bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet.get(), gcode.get(), sites.get(), params));
        //unique_ptr<bpp::FrequenciesSet> fSet(new bpp::FixedFrequenciesSet(model->getAlphabet(), model->getFrequencies()));
        //model_set.reset(bpp::SubstitutionModelSetTools::createHomogeneousModelSet(model.release(), fSet.release(), &tree));
        bpp_tree_like.reset(new bpp::RHomogeneousTreeLikelihood(tree, *sites, model.get(), rate_dist.get(), false, false, false));
    } else if(nh_opt == "general") {
        model_set.reset(bpp::PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet.get(), gcode.get(), sites.get(), params));
        bpp_tree_like.reset(new bpp::RNonHomogeneousTreeLikelihood(tree, *sites, model_set.get(), rate_dist.get(), false, false, false));
    } else throw std::runtime_error("Unknown non-homogeneous option: " + nh_opt);


    // lcfit-specific
    // Output
    string csv_like_path = bpp::ApplicationTools::getAFilePath("lcfit.output.likelihoods_file", params, true, false);
    ofstream csv_like_out(csv_like_path);
    string csv_ml_path = bpp::ApplicationTools::getAFilePath("lcfit.output.maxima_file", params, true, false);
    ofstream csv_ml_out(csv_ml_path);
    string csv_fit_path = bpp::ApplicationTools::getAFilePath("lcfit.output.fit_file", params, true, false);
    ofstream csv_fit_out(csv_fit_path);
    csv_fit_out << "node_id,branch_length,ll" << endl;;

    string csv_mltol_path = bpp::ApplicationTools::getAFilePath("lcfit.output.mltol_file", params, true, false);
    ofstream csv_mltol_out(csv_mltol_path);
    csv_mltol_out << "node_id,tolerance,ml_t,ml_est,n_eval,ml_brent,n_eval_brent" << endl;


    const vector<double> sample_points = bpp::ApplicationTools::getVectorParameter<double>(
            "lcfit.sample.branch.lengths",
            params, ',', "0.1,0.15,0.5,1.0");
    const vector<double> start = bpp::ApplicationTools::getVectorParameter<double>("lcfit.starting.values",
                                 params, ',', "1100,800,2.0,0.5");
    const vector<double> tolerance_values = bpp::ApplicationTools::getVectorParameter<double>(
            "lcfit.tolerance.values",
            params, ',', "1e-2,1e-3,1e-4,1e-5,1e-6");


    /*
     * Run evaluations on each node, write output
     */
    csv_like_out << "node_id,branch_length,bpp_ll,fit_ll" << endl;
    csv_ml_out << "node_id,t,t_hat,c,m,r,b" << endl;

    for (const int & node_id : tree.getNodesId()) {
        clog << "[lcfit eval] Node " << setw(4) << node_id << "\n";

        // Calculators
        TreeLikelihoodCalculator likelihood_calc(*bpp_tree_like);
        LCFitter fitter(start, sample_points, &likelihood_calc, &csv_fit_out);

        if (!tree.hasDistanceToFather(node_id)) 
	    continue;

        const bsm_t m = fitter.fit_model(node_id);
        const vector<LogLikelihoodComparison> evals = compare_log_likelihoods(tree, &likelihood_calc, node_id, m);

        // Write to CSV
        std::for_each(begin(evals), end(evals),
            [&csv_like_out](const LogLikelihoodComparison & e) { to_csv(csv_like_out, e); });

        const ModelFit fit = compare_ml_values(tree, node_id, m);
        to_csv(csv_ml_out, fit);

        double ml_est, ml_brent;
        size_t n_evals, n_evals_brent;
        for(const double tol : tolerance_values) {
            std::tie<double, size_t>(ml_est, n_evals) = fitter.estimate_ml_branch_length(node_id, tol);
            std::tie<double, size_t>(ml_brent, n_evals_brent) = fitter.estimate_ml_branch_length_brent(node_id, tol);
            csv_mltol_out << node_id
                << ',' << tol
                << ',' << tree.getDistanceToFather(node_id)
                << ',' << ml_est
                << ',' << n_evals
                << ',' << ml_brent
                << ',' << n_evals_brent
                << '\n';
        }
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
