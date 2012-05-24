#ifndef WEIGHTING_FUNCTIONS_H
#define WEIGHTING_FUNCTIONS_H

#include <climits>
#include <cassert>
#include <cmath>
//--
#include "CTF.h"

/**
 *  Weight functions for 8 bit pixel values.
 *  Includes the hat function in Debevec and Malik 1997, plus 
 *  additional functions.  Includes the ability to sample a weighting function
 *  into a LUT for fast evaluations.
 */
namespace WeightingFunctions{

    //Hat function
    //with default parameters, becomes the hat function in the debevec and Malik paper
    CTF::ctf_t hat(unsigned char value, unsigned char lower = 0, unsigned char upper = 255);
    void makeLUTHat(CTF::ctf_t* lut,  unsigned char lower = 0, unsigned char upper = 255);
}


//All weight functions are inlined---------------------------------------------

inline CTF::ctf_t WeightingFunctions::hat(unsigned char value, unsigned char lower, unsigned char upper){
    assert(pow(2,sizeof(unsigned char) * CHAR_BIT) == 256); //sanity check for 8-bit uchars
    assert(lower < upper); //Make sure hat bounds are correct

    const CTF::ctf_t z = static_cast<CTF::ctf_t>(value);
    const CTF::ctf_t Z_MIN = static_cast<CTF::ctf_t>(lower);
    const CTF::ctf_t Z_MAX = static_cast<CTF::ctf_t>(upper);
    return 
        (z <= Z_MIN || z >= Z_MAX) ? 0 : 
        ((z <= (Z_MIN + Z_MAX)/static_cast<CTF::ctf_t>(2.0)) ? (z - Z_MIN) : 
        (Z_MAX - z));
}





#endif //WEIGHTING_FUNCTIONS_H
