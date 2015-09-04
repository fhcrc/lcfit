#include <fstream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousTreeLikelihood.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

#include "lcfit.h"
#include "lcfit_cpp.h"
#include "lcfit_select.h"

typedef bpp::TreeTemplate<bpp::Node> Tree;

// Actual comparison code

/// Calculates the log-likelihood of a tree given data and a model
class TreeLikelihoodCalculator
{
  public:
    TreeLikelihoodCalculator(const bpp::TreeLikelihood& tl) :
        like(tl.clone())
    {
        like->initialize();
    }

    /// Calculate log-likelihood
    double calculate_log_likelihood()
    {
        return like->getLogLikelihood();
    }

    void set_branch_length(const size_t node, double length)
    {
        like->setParameterValue("BrLen" + std::to_string(node), length);
    }

    void get_branch_length(const size_t node) const
    {
        like->getParameterValue("BrLen" + std::to_string(node));
    }

  private:
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
    { }

    TreeLikelihoodCalculator* calc;
    const int node_id;
    size_t n_evals;

    double operator()(double d)
    {
        n_evals++;
        calc->set_branch_length(node_id, d);
        return calc->calculate_log_likelihood();
    }
};

double node_likelihood_callback(double t, void *data)
{
    NodeLikelihoodCalculator* calc = static_cast<NodeLikelihoodCalculator*>(data);
    return (*calc)(t);
}

/// Runs lcfit, generating parameters for the BSM
class LCFitter
{
  public:
    LCFitter(const std::vector<double>& start,
             const std::vector<double>& sample_points,
             TreeLikelihoodCalculator* calc) :
        start(start),
        sample_points(sample_points),
        calc(calc)
    { }

    /** Run lcfit, returning coefficients of the model.
     *
     *  Given a set of starting points, attempts to add
     *  (branch_length, ll) points until the function is
     *  non-monotonic, then fits the BSM using these sampled points.
     */
    bsm_t fit_model(const int node_id)
    {
        /* GSL requires at least four points for minimization. */
        if (sample_points.size() < 4) {
            throw std::invalid_argument("sample_points.size() < 4");
        }

        bsm_t m = {start[0], start[1], start[2], start[3]};
        NodeLikelihoodCalculator ll(calc, node_id);

        lcfit::LCFitResult fit_result =
                lcfit::fit_bsm_log_likelihood(ll, m, sample_points, 8);

        return fit_result.model_fit;
    }

    std::pair<double, size_t> estimate_ml_branch_length(const int node_id,
                                                        const double tol=1e-5,
                                                        bsm_t* model_out=nullptr)
    {
        /* GSL requires at least four points for minimization. */
        if (sample_points.size() < 4) {
            throw std::invalid_argument("sample_points.size() < 4");
        }

        bsm_t m = {start[0], start[1], start[2], start[3]};
        NodeLikelihoodCalculator ll(calc, node_id);

        log_like_function_t fn {node_likelihood_callback, &ll};
        bool success = false;

        std::vector<double> pts(sample_points);

        const double result = estimate_ml_t(&fn, pts.data(), pts.size(), tol,
                                            &m, &success, 1e-6, 1e4);

        // if (!success) {
        //     throw std::runtime_error("estimate_ml_t returned false");
        // }

        if (model_out) {
            *model_out = m;
        }

        return std::make_pair(result, ll.n_evals);
    }

  private:
    const std::vector<double> start;
    const std::vector<double> sample_points;
    TreeLikelihoodCalculator* calc;
};

int run_main(int argc, char** argv)
{
    bpp::BppApplication lcfit_bakeoff(argc, argv, "lcfit-bakeoff");

    std::map<string, string> params = lcfit_bakeoff.getParams();

    /********************/
    /* ARGUMENT PARSING */
    /********************/

    // Alphabet
    const std::unique_ptr<bpp::Alphabet> alphabet(bpp::SequenceApplicationTools::getAlphabet(params, "", false));

    // Genetic code
    std::unique_ptr<bpp::GeneticCode> gcode;
    if (bpp::CodonAlphabet* codon_alph = dynamic_cast<bpp::CodonAlphabet*>(alphabet.get())) {
        const string code_desc = bpp::ApplicationTools::getStringParameter("genetic_code", params, "Standard", "", true, true);
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
    } else if (nh_opt == "general") {
        model_set.reset(bpp::PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet.get(), gcode.get(), sites.get(), params));
        bpp_tree_like.reset(new bpp::RNonHomogeneousTreeLikelihood(tree, *sites, model_set.get(), rate_dist.get(), false, false, false));
    } else throw std::runtime_error("Unknown non-homogeneous option: " + nh_opt);


    // lcfit-specific
    string csv_fit_path = bpp::ApplicationTools::getAFilePath("lcfit.output.fit_file", params, true, false);
    ofstream csv_fit_out(csv_fit_path);

    const std::vector<double> sample_points = bpp::ApplicationTools::getVectorParameter<double>(
        "lcfit.sample.branch.lengths",
        params, ',', "0.1,0.15,0.5,1.0");
    const std::vector<double> start = bpp::ApplicationTools::getVectorParameter<double>("lcfit.starting.values",
                                                                                        params, ',', "1100,800,2.0,0.5");
    const std::vector<double> tolerance_values = bpp::ApplicationTools::getVectorParameter<double>(
        "lcfit.tolerance.values",
        params, ',', "1e-2,1e-3,1e-4,1e-5,1e-6");

    //
    // Run evaluations on each node, write output
    //

    csv_fit_out << "node_id,tolerance,n_evals,c,m,r,b,t_hat\n";
    csv_fit_out << std::setprecision(std::numeric_limits<double>::max_digits10);

    for (const int& node_id : tree.getNodesId()) {
        clog << "[lcfit eval] Node " << setw(4) << node_id << "\n";

        // Calculators
        TreeLikelihoodCalculator likelihood_calc(*bpp_tree_like);
        LCFitter fitter(start, sample_points, &likelihood_calc);

        if (!tree.hasDistanceToFather(node_id))
            continue;

        for (const double tol : tolerance_values) {
            std::stringstream ss;
            ss << std::scientific << std::setprecision(0) << tol;
            std::string tol_str = ss.str();
            std::string eval_key = "iterative:" + tol_str;

            bsm_t model_out;
            double t_hat;
            size_t n_evals;

            std::tie<double, size_t>(t_hat, n_evals) =
                    fitter.estimate_ml_branch_length(node_id, tol, &model_out);

            csv_fit_out << node_id << ","
                        << tol_str << ","
                        << n_evals << ","
                        << model_out.c << "," << model_out.m << ","
                        << model_out.r << "," << model_out.b << ","
                        << t_hat << "\n";
        }
    }

    return 0;
}

int main(int argc, char** argv)
{
    try {
        run_main(argc, argv);
    } catch (exception& e) {
        std::cerr << "Error:" << e.what() << std::endl;
        return 1;
    }
}
