#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>

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

const bpp::DNA DNA;
const bpp::RNA RNA;
const bpp::ProteicAlphabet AA;

template<typename T>
void print_vector(std::vector<T> v, const char delim='\t', ostream& out=std::cout) {
    for(auto &i : v) {
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
    std::unique_ptr<bpp::ISequence> reader = std::unique_ptr<bpp::ISequence>(
                fac.createReader(bpp::IOSequenceFactory::FASTA_FORMAT));
    std::unique_ptr<bpp::SequenceContainer> seqs = std::unique_ptr<bpp::SequenceContainer>(reader->read(in, alphabet));

    // Have to look up by name
    std::vector<std::string> names = seqs->getSequencesNames();
    bpp::SiteContainer *sequences = new bpp::VectorSiteContainer(alphabet);

    for(auto &name : names) {
        sequences->addSequence(seqs->getSequence(name), true);
    }

    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sequences);

    return sequences;
}

double get_ll(bpp::RHomogeneousTreeLikelihood like) {
    like.initialize();
    like.computeTreeLikelihood();
    return like.getLogLikelihood();
}

bpp::TreeTemplate<bpp::Node>* tree_of_newick_path(const std::string& path)
{
    bpp::Newick newick;
    return newick.read(path);
}

struct Options
{
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

int main(const int argc, const char **argv)
{
    Options options;
    try {
        options = parse_command_line(argc, argv);
    } catch(TCLAP::ArgException &e) {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }

    vector<double> t = {0.1, 0.2, 0.5, 1.}; // Branch lengths at which to sample.
    vector<double> l;
    vector<double> x = {1500, 1000, 2.0, 0.5}; // These are the starting values.

    // Reading the tree.
    unique_ptr<bpp::TreeTemplate<bpp::Node>> tree(tree_of_newick_path(options.newick_path));

    bpp::Newick newick;

    // Making the model.
    unique_ptr<bpp::SubstitutionModel> model(new bpp::JCnuc(&DNA));
    bpp::ConstantDistribution rate_dist(1.0, true);

    // Reading in alignment.
    ifstream in(options.alignment_path);
    unique_ptr<bpp::SiteContainer> input_aln (read_alignment(in, model->getAlphabet()));
    if(in.bad()) {
        cerr << "Cannot read from " << options.alignment_path << endl;
        return 1;
    }

    auto sons = tree->getRootNode()->getSonsId();
    //int to_change = sons.back();
    int to_change = 3;

    // Computing the tree likelihoods to be fit.
    for(int i = 0; i < t.size(); i++) {
        tree->setDistanceToFather(to_change, t[i]);
        //newick.write(*tree, cout);
        bpp::RHomogeneousTreeLikelihood like(*tree, *input_aln, model.get(), &rate_dist, false, false, false);
        l.push_back(get_ll(like));
    }

    int status = fit_ll(t.size(), t.data(), l.data(), x.data());

    cout << "fit values: ";
    print_vector(x);

    double c = x[0];
    double m = x[1];
    double r = x[2];
    double b = x[3];
    double t_hat = ml_t(c, m, r, b);

    printf("t_hat: %g\n", t_hat);

    ofstream file;

    file.open (options.output_path);
    double delta = 0.01;

    file << "#m=0,S=2" << endl;
    for(double t = delta; t <= 1.; t += delta) {
        tree->setDistanceToFather(to_change, t);
        //newick.write(*tree, cout);
        bpp::RHomogeneousTreeLikelihood like(*tree, *input_aln, model.get(), &rate_dist, false, false, false);
        file << t << " " << get_ll(like) << endl;
    }

    file << "#m=1,S=0" << endl;
    for(double t = delta; t <= 1.; t += delta) {
        file << t << " " << ll(t, c, m, r, b) << endl;
    }

    // The fit of the JC model.
    vector<double> lfit;
    for(int i = 0; i < t.size(); i++) {
        lfit.push_back(ll(t[i], c, m, r, b));
    }

    cout << "fit evaluations: " << endl;
    print_vector(l);
    print_vector(lfit);

    file.close();

    return status;
}
