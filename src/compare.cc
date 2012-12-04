#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <stdexcept>
#include <vector>

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/Io/IoSequenceFactory.h>
#include <Bpp/Seq/SiteTools.h>

#include "lcfit.h"

using namespace std;
typedef bpp::TreeTemplate<bpp::Node> Tree;

template<typename T>
void print_vector(std::vector<T> v, const char delim = '\t', ostream& out = std::cout)
{
    for(auto & i : v) {
        out << i << delim;
    }
    out << endl;
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

void to_csv(ostream& out, const Evaluation& e)
{
    out << e.node_id << ","
        << e.branch_length << ","
        << e.ll << ","
        << e.pred_ll << endl;
}

void to_csv(ostream& out, const ModelFit& f) {
    out << f.node_id << ","
        << f.t << ","
        << f.t_hat << ","
        << f.c << ","
        << f.m << ","
        << f.r << ","
        << f.b << endl;
}

int run_main(int argc, char** argv)
{
    bpp::BppApplication lcfit_compare(argc, argv, "lcfit-compare");

    map<string,string> params = lcfit_compare.getParams();

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

    // lcfit-specific
    // Output
    string csv_like_path = bpp::ApplicationTools::getAFilePath("lcfit.output.likelihoods_file", params, true, false);
    ofstream csv_like_out(csv_like_path);
    string csv_ml_path = bpp::ApplicationTools::getAFilePath("lcfit.output.maxima_file", params, true, false);
    ofstream csv_ml_out(csv_ml_path);
    string csv_fit_path = bpp::ApplicationTools::getAFilePath("lcfit.output.fit_file", params, true, false);
    ofstream csv_fit_out(csv_fit_path);
    csv_fit_out << "node_id,branch_length,ll" << endl;;

    const vector<double> sample_points = bpp::ApplicationTools::getVectorParameter<double>(
            "lcfit.sample.branch.lengths",
            params, ',', "0.1,0.15,0.5");
    const vector<double> start = bpp::ApplicationTools::getVectorParameter<double>("lcfit.starting.values",
            params, ',', "1500,1000,2.0,0.5");

    /*
     * Functions to run lcfit, evaluate results
     */
    // Calculate the log-likelihood of the tree in its current state
    auto get_ll = [&tree, &sites, &model, &rate_dist]() -> double {
        bpp::RHomogeneousTreeLikelihood like(tree, *sites, model.get(), rate_dist.get(), false, false, false);
        like.initialize();
        like.computeTreeLikelihood();
        return like.getLogLikelihood();
    };

    // Run lcfit, return coefficients of the model
    auto fit_model = [&tree,&start,&sample_points,&get_ll,&csv_fit_out](const int node_id) -> vector<double> {
        vector<double> l;
        vector<double> x = start;
        vector<double> t = sample_points;
        assert(is_sorted(sample_points.begin(), sample_points.end()));
        l.reserve(sample_points.size()+1);

        // Try starting points
        const double original_dist = tree.getDistanceToFather(node_id);
        for(const double &d : t) {
            tree.setDistanceToFather(node_id, d);
            l.push_back(get_ll());
        }

        // Add a fourth point based on the three sampled so far
        size_t offset; // Position
        double d;      // Branch length

        // Monotonically decreasing - add a point to the left
        if(l[0] > l[1] && l[1] > l[2]) {
            d = t[0] / 10.0; // TODO: What to choose?
            offset = 0;
        } else if (l[0] < l[1] && l[1] > l[2]) {
            // Maxima enclosed by (l[0], l[1])
            d = (t[1] + t[2]) / 2.0; // TODO: What to choose?
            offset = 2;
        } else {
            // Monotonically increasing
            d = (t[1] + t[2]) / 2.0; // TODO: What to choose?
            offset = 2;
        }
        t.insert(t.begin() + offset, d);
        tree.setDistanceToFather(node_id, d);
        l.insert(l.begin() + offset, get_ll());

        for(int i = 0; i < t.size(); ++i) {
            csv_fit_out << node_id << "," << t[i] << "," << l[i] << endl;
        }

        assert(is_sorted(sample_points.begin(), sample_points.end()));

        tree.setDistanceToFather(node_id, original_dist);
        const int status = fit_ll(t.size(), t.data(), l.data(), x.data());
        if(status) throw runtime_error("fit_ll returned: " + std::to_string(status));
        return x;
    };

    // Evaluate log-likelihood obtained by model fit versus actual log-likelihood at a variety of points
    auto evaluate_fit = [&tree,&get_ll](const int node_id, const vector<double>& x, const double delta) -> vector<Evaluation> {
        vector<Evaluation> evaluations;
        const double c = x[0];
        const double m = x[1];
        const double r = x[2];
        const double b = x[3];
        const double t = tree.getDistanceToFather(node_id);
        for(double t = delta; t <= 1.; t += delta) {
            tree.setDistanceToFather(node_id, t);
            double actual_ll = get_ll();
            double fit_ll = ll(t, c, m, r, b);
            evaluations.emplace_back(node_id, t, actual_ll, fit_ll);
        }
        tree.setDistanceToFather(node_id, t);
        return evaluations;
    };

    // Compare ML branch length estimate from lcfit to original branch length (presumed to be ML value)
    auto compare_ml_values = [&tree](const int node_id, const vector<double>& x) -> ModelFit {
        const double c = x[0];
        const double m = x[1];
        const double r = x[2];
        const double b = x[3];
        const double t_hat = ml_t(c, m, r, b);
        const double t = tree.getDistanceToFather(node_id);
        return ModelFit(node_id, t, t_hat, c, m, r, b);
    };

    /*
     * Run evaluations on each node, write output
     */
    csv_like_out << "node_id,branch_length,bpp_ll,fit_ll" << endl;
    csv_ml_out << "node_id,t,t_hat,c,m,r,b" << endl;
    for(const int& node_id : tree.getNodesId()) {
        cerr << "Node " << node_id << "\r";
        if(!tree.hasDistanceToFather(node_id)) continue;
        vector<double> x = fit_model(node_id);
        const vector<Evaluation> evals = evaluate_fit(node_id, x, 0.01);
        // Write to CSV
        std::for_each(begin(evals), end(evals),
                [&csv_like_out](const Evaluation& e) { to_csv(csv_like_out, e); });
        const ModelFit fit = compare_ml_values(node_id, x);
        to_csv(csv_ml_out, fit);
    }

    return 0;
}

int main(int argc, char **argv)
{
    try {
        run_main(argc, argv);
    } catch (exception& e) {
        cerr << "Error:" << e.what() << endl;
        return 1;
    }
}
