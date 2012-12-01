#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <vector>

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model/GTR.h>
#include <Bpp/Phyl/Model/HKY85.h>
#include <Bpp/Phyl/Model/JCnuc.h>
#include <Bpp/Phyl/Model/JTT92.h>
#include <Bpp/Phyl/Model/TN93.h>
#include <Bpp/Phyl/Model/WAG01.h>
#include <Bpp/Phyl/Model/AbstractSubstitutionModel.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Alphabet/RNA.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/Io/IoSequenceFactory.h>
#include <Bpp/Seq/SiteTools.h>

#include "tclap/CmdLine.h"

#include "lcfit.h"

using namespace std;
using Tree = bpp::TreeTemplate<bpp::Node>;

const bpp::DNA DNA;
const bpp::RNA RNA;
const bpp::ProteicAlphabet AA;

template<typename T>
void print_vector(std::vector<T> v, const char delim = '\t', ostream& out = std::cout)
{
    for(auto & i : v) {
        out << i << delim;
    }
    out << endl;
}

/// Read an alignment from a stream
/// \param in Input stream
/// \param alphabet The alphabet to use.
bpp::SiteContainer* read_alignment(std::istream &in, const bpp::Alphabet *alphabet)
{
    // Holy boilerplate - Bio++ won't allow reading FASTA files as alignments
    bpp::IOSequenceFactory fac;
    std::unique_ptr<bpp::ISequence> reader(
        fac.createReader(bpp::IOSequenceFactory::FASTA_FORMAT));
    std::unique_ptr<bpp::SequenceContainer> seqs(reader->read(in, alphabet));

    // Have to look up by name
    std::vector<std::string> names = seqs->getSequencesNames();
    bpp::SiteContainer *sequences = new bpp::VectorSiteContainer(alphabet);

    for(auto & name : names) {
        sequences->addSequence(seqs->getSequence(name), true);
    }

    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sequences);

    return sequences;
}


Tree* tree_of_newick_path(const std::string& path)
{
    bpp::Newick newick;
    return newick.read(path);
}

/// Command-line flags
struct Options {
    string newick_path;
    string alignment_path;
    string output_path;
};

Options parse_command_line(const int argc, const char **argv)
{
    TCLAP::CmdLine cmd("LCFIT", ' ', "probably not working");

    TCLAP::UnlabeledValueArg<string> newick_path(
        "newick_path", "Path to newick tree", true, "", "newick tree", cmd);
    TCLAP::UnlabeledValueArg<string> alignment_path(
        "alignment_path", "Path to alignment [fasta]", true, "", "fasta file", cmd);
    TCLAP::ValueArg<string> output_path(
        "o", "out", "Filename to write output values", false, "data.dat", "data file", cmd);

    cmd.parse(argc, argv);

    Options options;
    options.newick_path = newick_path.getValue();
    options.alignment_path = alignment_path.getValue();
    options.output_path = output_path.getValue();

    return options;
}

struct Evaluation {
    Evaluation(int node_id, double branch_length, double ll, double pred_ll) :
        node_id(node_id),
        branch_length(branch_length),
        ll(ll),
        pred_ll(pred_ll) {};
    int node_id;
    double branch_length, ll, pred_ll;
};

void to_csv(ostream& out, const Evaluation &e)
{
    out << e.node_id << ","
        << e.branch_length << ","
        << e.ll << ","
        << e.pred_ll << endl;
}

int main(const int argc, const char **argv)
{
    Options options;
    try {
        options = parse_command_line(argc, argv);
    } catch(TCLAP::ArgException &e) {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }

    const vector<double> t = {0.1, 0.2, 0.5, 1.}; // Branch lengths at which to sample.
    const vector<double> start = {1500, 1000, 2.0, 0.5}; // These are the starting values.

    // Reading the tree.
    unique_ptr<Tree> tree(tree_of_newick_path(options.newick_path));

    //bpp::Newick newick;

    // Making the model.
    unique_ptr<bpp::SubstitutionModel> model(new bpp::JCnuc(&DNA));
    bpp::ConstantDistribution rate_dist(1.0, true);

    // Reading in alignment.
    ifstream in(options.alignment_path);
    unique_ptr<bpp::SiteContainer> input_aln(read_alignment(in, model->getAlphabet()));
    if(in.bad()) {
        cerr << "Cannot read from " << options.alignment_path << endl;
        return 1;
    }

    // Calculate the log-likelihood of the tree in its current state
    auto get_ll = [&tree, &input_aln, &model, &rate_dist]() {
        bpp::RHomogeneousTreeLikelihood like(*tree, *input_aln, model.get(), &rate_dist, false, false, false);
        like.initialize();
        like.computeTreeLikelihood();
        return like.getLogLikelihood();
    };

    auto fit_model = [&](const int node_id) {
        vector<double> l;
        vector<double> x = start;
        vector<double> bls = t;
        l.reserve(x.size());
        const double original_dist = tree->getDistanceToFather(node_id);
        for(const double &d : t) {
            tree->setDistanceToFather(node_id, d);
            //newick.write(*tree, cout);
            l.push_back(get_ll());
        }
        tree->setDistanceToFather(node_id, original_dist);
        const int status = fit_ll(t.size(), bls.data(), l.data(), x.data());
        if(status) throw runtime_error("fit_ll returned: " + std::to_string(status));
        return x;
    };

    auto evaluate_fit = [&](const int node_id, const vector<double>& x, const double delta) {
        vector<Evaluation> evaluations;
        const double c = x[0];
        const double m = x[1];
        const double r = x[2];
        const double b = x[3];
        const double t_hat = ml_t(c, m, r, b);
        cerr << node_id << " t_hat: " << t_hat << endl;
        const double original_dist = tree->getDistanceToFather(node_id);
        for(double t = delta; t <= 1.; t += delta) {
            tree->setDistanceToFather(node_id, t);
            double actual_ll = get_ll();
            double fit_ll = ll(t, c, m, r, b);
            evaluations.emplace_back(node_id, t, actual_ll, fit_ll);
        }
        tree->setDistanceToFather(node_id, original_dist);
        return evaluations;
    };

    ofstream out(options.output_path);
    out << "node_id,branch_length,bpp_ll,fit_ll" << endl;
    for(const int& node_id : tree->getNodesId()) {
        if(!tree->hasDistanceToFather(node_id)) continue;
        vector<double> x = fit_model(node_id);
        vector<Evaluation> evals = evaluate_fit(node_id, x, 0.01);
        // Write to CSV
        std::for_each(begin(evals), end(evals),
                [&out](const Evaluation& e) { to_csv(out, e); });
    }

    return 0;
}
