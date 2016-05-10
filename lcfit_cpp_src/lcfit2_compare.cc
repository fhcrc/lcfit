#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include <Bpp/App/BppApplication.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

#include "gsl.h"
#include "lcfit2.h"

struct log_likelihood_data {
    bpp::TreeLikelihood* tl;
    int node_id;
};

double log_likelihood_callback(double t, void* data)
{
    if (t < 1e-6 || t > 1e4) { throw std::invalid_argument("t is out of bounds"); }

    log_likelihood_data* lnl_data = static_cast<log_likelihood_data*>(data);

    bpp::TreeLikelihood* tl = lnl_data->tl;
    int node_id = lnl_data->node_id;

    tl->setParameterValue("BrLen" + std::to_string(node_id), t);
    return tl->getLogLikelihood();
}

void compute_sampling_bounds(double (*lnl_fn)(double, void*), void* lnl_fn_args,
                             const double min_t, const double max_t, const double t0,
                             const double threshold, double* left_t, double* right_t)
{
    auto f = [lnl_fn, lnl_fn_args, threshold](double t) {
        return lnl_fn(t, lnl_fn_args) - threshold;
    };

    if (f(min_t) >= 0.0) {
        *left_t = min_t;
    } else {
        *left_t = gsl::find_root(f, min_t, t0);
    }

    if (f(max_t) >= 0.0) {
        *right_t = max_t;
    } else {
        *right_t = gsl::find_root(f, t0, max_t);
    }
}

void sample_curves(double (*lnl_fn)(double, void*), void* lnl_fn_args, const lcfit2_bsm_t* model,
                   const double min_t, const double max_t, const double t0,
                   const int node_id, std::ostream& output)
{
    const double lnl_t0 = lnl_fn(t0, lnl_fn_args);
    const double lnl_threshold = lnl_t0 - std::abs(0.01 * lnl_t0);

    double left_t;
    double right_t;

    compute_sampling_bounds(lnl_fn, lnl_fn_args, min_t, max_t, t0,
                            lnl_threshold, &left_t, &right_t);

    std::cerr << "left = " << left_t << ", right = " << right_t << "\n";

    const size_t n_samples = 501;
    const double delta = (right_t - left_t) / (n_samples - 1);

    for (size_t i = 0; i < n_samples; ++i) {
        const double t = left_t + (i * delta);
        const double empirical_lnl = lnl_fn(t, lnl_fn_args);
        const double fit_lnl = lcfit2_lnl(t, model);

        output << node_id << ","
               << t << ","
               << empirical_lnl << ","
               << fit_lnl << "\n";
    }
}

int run_main(int argc, char** argv)
{
    bpp::BppApplication lcfit2_compare(argc, argv, "lcfit2-compare");

    std::map<std::string, std::string> params = lcfit2_compare.getParams();

    //
    // Argument parsing
    //

    // Alphabet
    const std::unique_ptr<bpp::Alphabet> alphabet(bpp::SequenceApplicationTools::getAlphabet(params, "", false));

    // Genetic code
    std::unique_ptr<bpp::GeneticCode> gcode;
    if (bpp::CodonAlphabet* codon_alph = dynamic_cast<bpp::CodonAlphabet*>(alphabet.get())) {
        const std::string code_desc = bpp::ApplicationTools::getStringParameter("genetic_code", params, "Standard", "", true, true);
        bpp::ApplicationTools::displayResult("Genetic Code", code_desc);
        gcode.reset(bpp::SequenceApplicationTools::getGeneticCode(codon_alph->getNucleicAlphabet(), code_desc));
    }

    // Sites
    std::unique_ptr<bpp::VectorSiteContainer> all_sites(bpp::SequenceApplicationTools::getSiteContainer(alphabet.get(), params));
    std::unique_ptr<bpp::VectorSiteContainer> sites(bpp::SequenceApplicationTools::getSitesToAnalyse(*all_sites, params, "", true, false));
    all_sites.reset();
    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    // Tree
    std::unique_ptr<bpp::Tree> in_tree(bpp::PhylogeneticsApplicationTools::getTree(params));
    bpp::TreeTemplate<bpp::Node> tree(*in_tree);

    // Rate dist
    std::unique_ptr<bpp::DiscreteDistribution> rate_dist(bpp::PhylogeneticsApplicationTools::getRateDistribution(params));

    // Model
    // See bppSeqGen.cpp L253
    const std::string nh_opt = bpp::ApplicationTools::getStringParameter("nonhomogeneous", params, "no", "", true, false);
    std::unique_ptr<bpp::SubstitutionModel> model;
    std::unique_ptr<bpp::SubstitutionModelSet> model_set;
    std::unique_ptr<bpp::TreeLikelihood> bpp_tree_like;

    if (nh_opt == "no") {
        model.reset(bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet.get(), gcode.get(), sites.get(), params));
        //std::unique_ptr<bpp::FrequenciesSet> fSet(new bpp::FixedFrequenciesSet(model->getAlphabet(), model->getFrequencies()));
        //model_set.reset(bpp::SubstitutionModelSetTools::createHomogeneousModelSet(model.release(), fSet.release(), &tree));
        bpp_tree_like.reset(new bpp::RHomogeneousTreeLikelihood(tree, *sites, model.get(), rate_dist.get(), false, false, false));
    }
    else if (nh_opt == "general") {
        model_set.reset(bpp::PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet.get(), gcode.get(), sites.get(), params));
        bpp_tree_like.reset(new bpp::RNonHomogeneousTreeLikelihood(tree, *sites, model_set.get(), rate_dist.get(), false, false, false));
    }
    else {
        throw std::runtime_error("Unknown non-homogeneous option: " + nh_opt);
    }

    // Output files
    std::string lnl_filename = bpp::ApplicationTools::getAFilePath("lcfit2.output.lnl_file", params, true, false);
    std::ofstream lnl_output(lnl_filename);
    lnl_output << "node_id,t,empirical,lcfit2\n";
    lnl_output << std::setprecision(std::numeric_limits<double>::max_digits10);

    std::string lcfit2_filename = bpp::ApplicationTools::getAFilePath("lcfit2.output.fit_file", params, true, false);
    std::ofstream lcfit2_output(lcfit2_filename);
    lcfit2_output << "node_id,c,m,t0,d1,d2\n";
    lcfit2_output << std::setprecision(std::numeric_limits<double>::max_digits10);

    //
    // Run evaluations on each node, write output
    //

    for (const int& node_id : tree.getNodesId()) {
        if (!tree.hasDistanceToFather(node_id)) {
            continue;
        }

        clog << "[lcfit2 eval] Node " << std::setw(4) << node_id << "\n";

        //
        // initialize the tree likelihood calculator
        //

        // clone the tree likelihood calculator to ensure that any
        // changes are isolated to computations for this node

        std::unique_ptr<bpp::TreeLikelihood> tl(bpp_tree_like->clone());

        // derivatives can be enabled/disabled before the call to
        // `tl->initialize()`. if they are disabled, computed
        // derivative values will be zero. it appears that "enabled"
        // is the default, but we enable them explicitly for clarity.

        tl->enableDerivatives(true);
        tl->initialize();

        log_likelihood_data lnl_data = { tl.get(), node_id };

        //
        // find t0
        //

        bpp::DiscreteRatesAcrossSitesTreeLikelihood* drastl = dynamic_cast<bpp::DiscreteRatesAcrossSitesTreeLikelihood*>(tl.get());

        bpp::ParameterList pl;
        pl.addParameter(drastl->getBranchLengthsParameters().getParameter("BrLen" + std::to_string(node_id)));
        bpp::OptimizationTools::optimizeBranchLengthsParameters(drastl, pl);

        const double t0 = tl->getParameterValue("BrLen" + std::to_string(node_id));

        std::cerr << "t0 = " << t0 << "\n";

        //
        // find d1 and d2 at t0
        //

        // ensure branch length is set to t0
        tl->setParameterValue("BrLen" + std::to_string(node_id), t0);

        // GOTCHA: for unknown reasons the values returned by these
        // functions appear to be the negatives of the corresponding
        // derivatives, so we negate them again to get the proper
        // values.
        const double d1 = -(tl->getFirstOrderDerivative("BrLen" + std::to_string(node_id)));
        const double d2 = -(tl->getSecondOrderDerivative("BrLen" + std::to_string(node_id)));

        std::cerr << "d1(t0) = " << d1 << "\n";
        std::cerr << "d2(t0) = " << d2 << "\n";

        //
        // fit with lcfit2
        //

        const double min_t = 1e-6;
        const double max_t = 1e4;

        lcfit2_bsm_t model = {1100.0, 800.0, t0, d1, d2};

        if (std::abs(d1) < 0.1) {
            const double alpha = 0.0;

            // GOTCHA: this function will change the current branch length
            lcfit2_fit_auto(&log_likelihood_callback, &lnl_data, &model, min_t, max_t, alpha);

            lcfit2_output << node_id << "," << model.c << "," << model.m << ","
                          << model.t0 << "," << model.d1 << "," << model.d2 << "\n";
        }

        //
        // sample empirical and lcfit2 curves
        //

        if (std::abs(d1) < 0.1) {
            // GOTCHA: this function will change the current branch length
            sample_curves(&log_likelihood_callback, &lnl_data, &model,
                          min_t, max_t, t0, node_id, lnl_output);
        }
    }

    return 0;
}

int main(int argc, char** argv)
{
    try {
        run_main(argc, argv);
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
