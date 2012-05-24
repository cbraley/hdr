#include "WeightingFunctions.h"

void WeightingFunctions::makeLUTHat(CTF::ctf_t* lut,  unsigned char lower, unsigned char upper){
    assert(lut != NULL);
    for(int pixVal = 0; pixVal < 256; pixVal++){
        lut[pixVal] = WeightingFunctions::hat(static_cast<unsigned char>(pixVal), lower, upper);
    }
}

