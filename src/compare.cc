#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model/GTR.h>
#include <Bpp/Phyl/Model/HKY85.h>
#include <Bpp/Phyl/Model/JCnuc.h>
#include <Bpp/Phyl/Model/JTT92.h>
#include <Bpp/Phyl/Model/TN93.h>
#include <Bpp/Phyl/Model/WAG01.h>
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

#include "lcfit.h"

using namespace std;

const bpp::DNA DNA;
const bpp::RNA RNA;
const bpp::ProteicAlphabet AA;

template<typename T> void print_vector(std::vector<T> v) {
    for(auto it = v.begin(); it != v.end(); ++it) {
        std::cout << *it << "\t";
    }
    std::cout << endl;
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

    for(auto name : names) {
        sequences->addSequence(seqs->getSequence(name), true);
    }

    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sequences);

    return sequences;
}


int main(void)
{
    const size_t n = 3;
    double t[n] = {0.25, 0.5, 1}; // Branch lengths at which to sample.
    vector<double> l;
    double x[n] = { 100.0, 2.0, 1.0 }; // These are the starting values.
    std::string aln_fname = "sts/data/test.fasta";

    // Setting up output.
    //string output_filename = output_path.getValue();
    ostream *output_stream;
    //ofstream output_ofstream;
    //if(output_filename == "-") {
        output_stream = &cout;
    //} else {
    //    output_ofstream.open(output_filename);
    //    output_stream = &output_ofstream;
    //}


    // Reading the tree.
    std::string tree_str = "(F:0.18861,(D:0.19026,E:0.14671)0.999:0.29165,(H:0.29607,(G:0.13891,(C:0.15128,(A:0.21353,B:0.20443)0.998:0.29027)0.936:0.17926)1.000:0.50775)0.837:0.10589);";
    std::unique_ptr<bpp::TreeTemplate<bpp::Node>> tree(bpp::TreeTemplateTools::parenthesisToTree(tree_str));

    // Making the model.
    shared_ptr<bpp::SubstitutionModel> model;
    model = make_shared<bpp::JCnuc>(&DNA);
    bpp::ConstantDistribution rate_dist(1.0, true);

    // Reading in alignment.
    ifstream in(aln_fname);
    unique_ptr<bpp::SiteContainer> input_aln (read_alignment(in, model->getAlphabet()));
    if(in.bad()) {
        cerr << "Cannot read from " << aln_fname << endl;
        return 1;
    }

    // Initializing the likelihood computation engine.
    bpp::RHomogeneousTreeLikelihood like(*tree, *input_aln, model.get(), &rate_dist, false, false, false);
    like.initialize();

    // Computing the tree likelihoods to be fit.
    for(int i = 0; i < n; i++) {
        like.computeTreeLikelihood();
        l.push_back(like.getLogLikelihood());
    }

    print_vector(l);

    int status = fit_ll(n, t, l.data(), x);

    printf("fit values: ");
    for(int i = 0; i < n; i++)
        printf("%g\t", x[i]);
    printf("\n");

    double c = x[0];
    double r = x[1];
    double m = x[2];
    double t_hat = ml_t(c, m, r);

    printf("t_hat: %g\n", t_hat);

    return status;
}

