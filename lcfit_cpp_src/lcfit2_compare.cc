#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include <Bpp/App/BppApplication.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

#include "lcfit2.h"

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

    //
    // Run evaluations on each node, write output
    //

    for (const int& node_id : tree.getNodesId()) {

        // Calculators
        //TreeLikelihoodCalculator likelihood_calc(*bpp_tree_like);
        //LCFitter fitter(start, sample_points, &likelihood_calc, &csv_fit_out);

        // pseudocode equivalent of the initialization steps that
        // occur above to obtain a log-likelihood calculator. this
        // should happen somewhere else and must ensure that any
        // changes are isolated to computations for this node

        std::unique_ptr<bpp::TreeLikelihood> tl(bpp_tree_like->clone());
        tl->initialize();

        auto lnl_callback = [node_id](double t, void* data) {
            if (t < 1e-6 || t > 1e4) { throw std::invalid_argument("t is out of bounds"); }

            bpp::TreeLikelihood* calc = static_cast<bpp::TreeLikelihood*>(data);
            calc->setParameterValue("BrLen" + std::to_string(node_id), t);
            return calc->getLogLikelihood();
        };

        // should we/where do we enable derivative computation? can we
        // turn them on when we need them and disable them when we
        // don't, or is it a one-time thing?
        //
        // virtual void bpp::TreeLikelihood::enableDerivatives(bool yn)

        // end pseudocode

        clog << "[lcfit eval] Node " << std::setw(4) << node_id << "\n";

        if (!tree.hasDistanceToFather(node_id)) {
            continue;
        }

        //
        // sample empirical curve
        //

        // test: compute log-likelihood at the current branch length
        const double t = tl->getParameterValue("BrLen" + std::to_string(node_id));
        const double lnl = lnl_callback(t, tl.get());

        std::cerr << "ell(" << t << ") = " << lnl << "\n";

        // ...

        // find t0
        //
        // is the maximum-likelihood branch length just the value
        // that's already attached to the branch when the tree is
        // loaded? if so we don't actually need to do anything except
        // cache the value before we start sampling the curve

        // ...

        //
        // find d1 and d2
        //
        // it may be possible to use Bio++ methods to get these values:
        //
        // virtual double getFirstOrderDerivative (const std::string &variable)
        // virtual double getSecondOrderDerivative (const std::string &variable)


        // ...

        //
        // fit with lcfit2
        //

        // ...

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
