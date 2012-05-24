#ifndef CTF_SOLVER_H
#define CTF_SOLVER_H

#include <iostream>
#include <string>
#include <vector>
#include <cassert>
//--
#include "CTF.h"
#include "WeightingFunctions.h"

//TODO: Account for blooming pixels
class CTFSolver{
public:

    //Structure containing a path to an image along with the exposure time
    //for the image
    typedef struct ImageExposurePair{

        ImageExposurePair(long microSeconds, const std::string& path);

        long microseconds;     //Exposure time in microseconds
        std::string imagePath; //Path to LDR image on disk

        //Get floating point time value
        CTF::ctf_t getTime()const{ return (CTF::ctf_t) (microseconds); }

        //Comparison based on exposure time
        bool operator<(const ImageExposurePair& other)const;
        
        friend std::ostream& operator<<(std::ostream& os, const ImageExposurePair& im);
    }ImageExposurePair;

    /**
     *  Construct a CTF solver.  The actual image data on disk is lazy-loaded 
     *  as needed in an effort to conserve memory.  The loading occurs once 
     *  solveCTF() is called.
     *
     *  @param images is a list of ImageExposurePair structs.  Each image, as specified 
     *  by the imagePath member in struct ImageExposurePair, MUST exist on disk.  This is
     *  not checked by the code, and code will fail on solveCTF() if this is not the case.
     *  Furthermore, all images must have the same dimensions.
     *  @param numSamps is the number of random image samples to take.  Defaults to 1000.
     *   Should always be greater than 256 in order to make the resulting linear system
     *   at least square.  Values > 256 result in an overdetermined system.
     *  @param smoothingParam corresonds to lambda in the Debevec and Malik paper.  Defaults to 1.0.
     *  @param channel is which color channel of the image we should use.  Defaults to 0, thus being 
     *  appropriate for monochrome images by default.  To use the green channel of a RGB camera, set
     *  channel = 1.
     */
    CTFSolver(const std::vector<ImageExposurePair>& images, 
        size_t numSamps = 1000,
        CTF::ctf_t smoothingParam = static_cast<CTF::ctf_t>(1.0),
        size_t channel = 0);


    typedef struct PixelResult{
        CTF::ctf_t irradiance;
        int x,y;
    }PixelResult;

    bool writePixelPoints(const std::vector<PixelResult>& pixels,
        std::ostream& os)const;

    CTF solve(std::vector<PixelResult>* retPixels = NULL)const; 

    friend std::ostream& operator<<(std::ostream& os, const CTFSolver& solver);

    //Types of possible weighting functions
    enum WeightingFunc{HAT, HAT_10};
    void setWeightingFunc(WeightingFunc func);
    WeightingFunc getWeightingFunc()const;

    void setNumImageSamples(size_t numSamps);
    size_t getNumImageSamples()const;

    void setSmoothingValue(CTF::ctf_t lambdaVal);
    CTF::ctf_t getSmoothingValue()const;

    size_t getChannelIndex()const;
    void setChannelIndex(size_t chanIndex);

    /// Load a particular HDR stack and make sure all the images exist.
    /// Also, return the width, height, and maximum color channels in the stack.
    ///
    /// @return true if all images exist and have the same dimensions.  Returns 
    /// false otherwise.
    /// @outWidth is set to the width of the images.
    /// @outHeight is set to the width of the images.
    /// @outMinNumChans is set to the minimum number of color channels found.
    /// @reason is a return parameter that describes errors(if any)
    static bool checkImagesOK(std::vector<CTFSolver::ImageExposurePair>& images,
        int& outWidth, int& outHeight, int& outMinNumChans, std::string* reason = NULL);


private:
    //Non-Copyable
    CTFSolver(const CTFSolver& other);
    CTFSolver& operator=(const CTFSolver& rhs);

    //Data
    std::vector<ImageExposurePair> imdata; //List of images along with exposure times
    CTF::ctf_t lambda; //Smoothing value
    size_t chan; //Color channel index
    size_t numSamples; //How many random samples to take from the image?
    WeightingFunc wFunc; //Which weighting function are we using?

    //Helper functions
    CTF::ctf_t hatFunc(unsigned char zVal)const;
    CTF::ctf_t hatFuncParameterized(unsigned char zVal, unsigned char cut)const;
};

inline void CTFSolver::setWeightingFunc(WeightingFunc func){
    wFunc = func;
}
inline CTFSolver::WeightingFunc CTFSolver::getWeightingFunc()const{
    return wFunc;
}

inline void CTFSolver::setNumImageSamples(size_t numSamps){
    numSamples = numSamps;
}
inline size_t CTFSolver::getNumImageSamples()const{
    return numSamples;
}

inline void CTFSolver::setSmoothingValue(CTF::ctf_t lambdaVal){
    lambda = lambdaVal;
}
inline CTF::ctf_t CTFSolver::getSmoothingValue()const{
    return lambda;
}

inline size_t CTFSolver::getChannelIndex()const{
    return chan;
}
inline void CTFSolver::setChannelIndex(size_t chanIndex){
    chan = chanIndex;
}

inline CTF::ctf_t CTFSolver::hatFunc(unsigned char zVal)const{
    /*
    const CTF::ctf_t z = static_cast<CTF::ctf_t>(zVal);
    const CTF::ctf_t Z_MIN = static_cast<CTF::ctf_t>(0.0  );
    const CTF::ctf_t Z_MAX = static_cast<CTF::ctf_t>(255.0);
    assert(z >= Z_MIN && z <= Z_MAX);
    return 
        ((z <= (Z_MIN + Z_MAX)/static_cast<CTF::ctf_t>(2.0)) ? (z - Z_MIN) : 
        (Z_MAX - z));
    */
    return WeightingFunctions::hat(zVal);
}


inline CTF::ctf_t CTFSolver::hatFuncParameterized(unsigned char zVal,
    unsigned char cut)const
{
    /*
    assert(cut < 255);

    const CTF::ctf_t z = static_cast<CTF::ctf_t>(zVal);
    const CTF::ctf_t Z_MIN = static_cast<CTF::ctf_t>(cut  );
    const CTF::ctf_t Z_MAX = static_cast<CTF::ctf_t>(255.0 - cut);
    return (z <= Z_MIN || z >= Z_MAX) ? static_cast<CTF::ctf_t>(0.0) : 
        (
        ((z <= (Z_MIN + Z_MAX)/static_cast<CTF::ctf_t>(2.0)) ? (z - Z_MIN) : 
        (Z_MAX - z))
        );
    */
    return WeightingFunctions::hat(zVal, cut, 255-cut);
}


inline bool CTFSolver::ImageExposurePair::operator<(const ImageExposurePair& other)const{
    return microseconds < other.microseconds;
}


#endif //CTF_SOLVER_H
