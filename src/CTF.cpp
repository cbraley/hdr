#include "CTF.h"
//--
#include <fstream>
#include <cstdlib>

CTF::CTF(){
    CTF::ctf_t v          = static_cast<CTF::ctf_t>(0.0);
    const CTF::ctf_t step = static_cast<CTF::ctf_t>(0.0/256.0);
    for(int i = 0; i < 256; i++){
        data.push_back(v);
        v += step;
    }

}

CTF::CTF(const std::vector<ctf_t>& values) : 
    data(values)
{
    assert(values.size() == 256);
}

std::ostream& operator<<(std::ostream& os, const CTF& ctf){
    for(size_t i = 0; i < ctf.data.size(); i++){
        os << ctf.data[i] << std::endl;
    }
    return os;
}


bool CTF::loadCTF(CTF& ctf, const std::string& fileName){
    std::fstream fs(fileName.c_str(), std::fstream::in);
    if(! fs.good() ){
        fs.close();
        return false;
    }

    //Read the data
    assert(ctf.data.size() == 256);
    std::string str;
    for(int i = 0; i < 256; i++){
        std::getline(fs, str);

        //Convert to ctf_t
        CTF::ctf_t exposure = static_cast<CTF::ctf_t>(strtod(str.c_str(), NULL));

        //exposure  == 0.0 indicates a parse error (strtod's way of error handling)
        //exposure < 0 indicates a physically implausable response curve
        if(exposure <= static_cast<CTF::ctf_t>(0.0) ){
            fs.close();
            return false;
        }
        
        //Update the CTF
        ctf.data[i] = exposure;
    }


    const bool worked = fs.good();
    fs.close();
    return worked;
}


CTF CTF::makeLinearCTF(CTF::ctf_t maxCTFValue,
    CTF::ctf_t minCTFValue)
{
    assert(minCTFValue < maxCTFValue);

    CTF ctf;
    CTF::ctf_t step = (maxCTFValue - minCTFValue) / static_cast<CTF::ctf_t>(255);
    CTF::ctf_t val = minCTFValue;
    for(int pix = 0; pix <= 255; pix++){
        ctf.data[pix] = val; 
        val += step;
    }

    return ctf;
}
