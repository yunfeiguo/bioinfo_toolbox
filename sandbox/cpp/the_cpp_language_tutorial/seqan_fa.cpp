//
// Created by guoy28 on 3/14/17.
//

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include "seqan_fa.h"
void seqan_fa_example() {
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::CharString seqFileName = "/Users/guoy28/Downloads/bioinfo_toolbox/sandbox/cpp/the_cpp_language_tutorial/resources/example.fa";

    try {
        seqan::SeqFileIn seqFileIn(toCString(seqFileName));
        readRecord(id, seq, seqFileIn);
        std::cout << id << '\t' << seq << '\n';
    } catch (std::exception const &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }

    seqan::StringSet <seqan::CharString> ids;
    seqan::StringSet <seqan::Dna5String> seqs;
    std::vector<std::string> mySeqs;
    std::vector<std::string> myIds;

    try {
        seqan::SeqFileIn seqFileIn1(seqan::toCString(seqFileName));
        seqan::readRecords(ids, seqs, seqFileIn1);
    } catch (std::exception const &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }

    for (unsigned i = 0; i < seqan::length(ids); i++) {
        std::cout << ids[i] << '\t' << seqs[i] << '\n';
        std::string s = to_string(seqs[i]);
        std::string k(seqan::toCString(ids[i]));
        mySeqs.push_back(s);
        myIds.push_back(k);
    }
}
