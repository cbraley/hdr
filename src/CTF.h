#ifndef CTF_H
#define CTF_H

#include <vector>
#include <string>
#include <iostream>
#include <cassert>

/**
 *  Class representing the camera transfer function of a camera(the CTF).  The goal of this
 *  software is to solve for these curves.  The CTF is a function, internal to the camera, that
 *  converts "irradiance" to 8 bit values.
 */
class CTF{
public:

    //Single precision(float) should typically be sufficient, but you
    //could change this to double if you like
    typedef float ctf_t;

    /**
     *  Initialize a CTF.
     */
    CTF(const std::vector<ctf_t>& values);

    CTF();

    /**
     *  Load the CTF value for a particular pixel value in the range [0, 255].
     */
    ctf_t operator()(unsigned char pixelVal)const;

    /**
     *  Write CTF to a stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const CTF& ctf);

    /**
     *  Load a CTF from disk.  Return true on success, false on failure.
     *  The CTF ctf is modified to return the result.
     *  This is not a constructor because constructors have no good way to indicate failure.
     */
    static bool loadCTF(CTF& ctf, const std::string& fileName);

    /**
     *  Create a linear CTF.  
     *
     *  @param maxCTFValue is the maximum CTF value.
     *  @param minimumCTFValue is the minimum CTF value.  Defaults to 0.
     */
    static CTF makeLinearCTF(CTF::ctf_t maxCTFValue,
        CTF::ctf_t minCTFValue = static_cast<CTF::ctf_t>(0.0));


private:

    std::vector<CTF::ctf_t> data; //Array of length num_bits (currently 256)
};


inline CTF::ctf_t CTF::operator()(unsigned char pixelVal)const{
    assert(data.size() == 256);
    return data[pixelVal];
}



#endif //CTF_H
