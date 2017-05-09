//
// Created by guoy28 on 3/14/17.
//

#ifndef THE_CPP_LANGUAGE_TUTORIAL_SEQAN_FA_H
#define THE_CPP_LANGUAGE_TUTORIAL_SEQAN_FA_H


class seqan_fa {

};

void seqan_fa_example();

inline std::string seqanBase2String(seqan::Dna5 base){

    //std::cout<<static_cast<unsigned char>(base)<<" :static cast: ";
    //std::cout<<(unsigned)seqan::ordValue(base)<<std::endl;

    switch((unsigned)seqan::ordValue(base))
    {
        case 0: return "A";
        case 1: return "C";
        case 2: return "G";
        case 3: return "T";
        default: return "N";
    }
    //should never reach here
    return "";

}

inline std::string to_string(const seqan::Dna5String& seq){
    std::string result;
    for (unsigned int i = 0; i < seqan::length(seq); ++i) {
        result += seqanBase2String(seq[i]);
    }
    return result;
}


#endif //THE_CPP_LANGUAGE_TUTORIAL_SEQAN_FA_H
