#include "CTFSolver.h"
//--
#include <cassert>
#include <cstdlib>
//--
#define cimg_display 0
#include "cimg/CImg.h"
#undef cimg_display
using namespace cimg_library;
//--
#include "Eigen/Eigen"
#include "Eigen/Dense"


CTFSolver::ImageExposurePair::ImageExposurePair(long microSeconds, const std::string& path) : 
    microseconds(microSeconds), imagePath(path)
{}


//TODO: Use a better RNG
static int randomInt(int minInclusive, int maxExclusive){
    assert(maxExclusive > minInclusive);
    const int val =  minInclusive + (rand() % (maxExclusive - minInclusive));
    assert(val >= minInclusive && val < maxExclusive);
    return val;
}


typedef struct SamplePos{
    SamplePos(int xPos = -1, int yPos = -1) : 
        x(xPos), y(yPos){}
    int x,y;
}SamplePos;

static std::vector<SamplePos> genRandomSamples(int widthMax, int heightMax, int numSamps){

    //Allocate space
    std::vector<SamplePos> ret;
    ret.resize((size_t)numSamps);
    
    //Add samples
    //TODO: Use stochastic sampling not pure random
    for(int i = 0; i < numSamps; i++){
        SamplePos sample(randomInt(0,widthMax), randomInt(0,heightMax));
        ret[i] = sample;
    }

    return ret;
}


CTFSolver::CTFSolver(const std::vector<ImageExposurePair>& images, 
    size_t numSamps,
    CTF::ctf_t smoothingParam,
    size_t channel) : 
    imdata(images), lambda(smoothingParam), chan(channel), numSamples(numSamps),
    wFunc(HAT)
{
    assert(images.size() >= 2);
    //assert(numSamples > 256);
}


bool CTFSolver::writePixelPoints(const std::vector<PixelResult>& pixels,
    std::ostream& os)const
{
    for(size_t j = 0; j < imdata.size(); j++){ //Loop over images

        //Load current image
        CImg<unsigned char> currIm(imdata[j].imagePath.c_str());
        
        //Loop over samples
        for(size_t i = 0; i < pixels.size(); i++){
            const int x = pixels[i].x;
            const int y = pixels[i].y;
            const CTF::ctf_t irradiance = pixels[i].irradiance;
            const unsigned char pixelValue = currIm(x,y,0,chan);
            const CTF::ctf_t exposure = imdata[j].getTime() * irradiance;
            os << static_cast<int>(pixelValue) << "     " << exposure << std::endl;
        }
    }


    return os.good();
}

CTF CTFSolver::solve(std::vector<PixelResult>* retPixels)const{

    //n = 256 for 8 bit images
    const int n = 256;

    //Make a lookup table for our weighting function
    //We don't care that we branch inside of the loop since this LUT is created once
    CTF::ctf_t wLut[n];
    for(int i = 0; i < n; i++){
        CTF::ctf_t val = static_cast<CTF::ctf_t>(1.0);
        const unsigned char pixVal = static_cast<unsigned char>(i);
        switch(wFunc){
            case HAT:
                val = hatFunc(pixVal);
                break;
            case HAT_10:
                val = hatFuncParameterized(pixVal, 10);
                break;
            default:
                //This case should never occur, since the switch statement 
                //should be exhaustive for all possible weighting functions
                assert(false);
                val = hatFunc(pixVal);
                break;
        }
        wLut[i] = val;
    }

    //Find dimensions of first image
    CImg<unsigned char> firstIm(imdata[0].imagePath.c_str());
    const int firstWidth  = firstIm.width();
    const int firstHeight = firstIm.height();

    //Generate vector of random sample positions
    std::vector<SamplePos> samplePositions = genRandomSamples(firstWidth, firstHeight, numSamples);
    assert(samplePositions.size() == (size_t)numSamples);

    //Allocate space for solving linear system

    //Create left hand side matrix A in Ax=b
    Eigen::Matrix<CTF::ctf_t, Eigen::Dynamic, Eigen::Dynamic> A(
        numSamples * imdata.size() + n + 1, //rows
        n + numSamples                      //columns
        );
    //Right hand side vector
    Eigen::Matrix<CTF::ctf_t, Eigen::Dynamic, Eigen::Dynamic> b(
        numSamples * imdata.size() + n + 1, //rows
        1
    );

    //Populate the linear system

    //   Fitting equations
    //   Note that the loop order here is reversed as compared to the Debevec paper to allow
    //   us to load 1 image at a time and save memory
    //   TODO: Catch loading errors
    int k = 0;
    for(size_t j = 0; j < imdata.size(); j++){ //Loop over images

        //Load current image
        CImg<unsigned char> currIm(imdata[j].imagePath.c_str());
        assert(currIm.width() == firstWidth);
        assert(currIm.height() == firstHeight);
        assert((size_t)currIm.spectrum() > chan);

        for(size_t i = 0; i < numSamples; i++){ //Loop over sample positions
            
            //Get pixel value at current sample
            const int x = samplePositions[i].x;
            const int y = samplePositions[i].y;
            const unsigned char pixVal = currIm(x,y,0,chan);

            //Get value of weighting function
            const CTF::ctf_t w = wLut[pixVal];
            
            //Update A matrix
            A(k,pixVal)  = w;
            A(k,n+i)     = -w;

            //Update RHS b vector
            const CTF::ctf_t t = imdata[j].getTime();
            assert(t > 0.0);

            b(k,0)        = w * log(t);

            ++k;
        }
    }

    //    Fix the curve
    A(k++, 128) = static_cast<CTF::ctf_t>(1.0);


    //   Include regularization
    for(int i = 0; i <= n-2; i++){
        const int lval = i+1;
        assert(lval >= 0 && lval <= 255);
        const CTF::ctf_t w = wLut[lval];
        A(k,i  ) = lambda * w;
        A(k,i+1) = -2.0 * lambda * w;
        A(k,i+2) = lambda * w;

        ++k;
    }


    //At long last, solve the system
    Eigen::Matrix<CTF::ctf_t, Eigen::Dynamic, Eigen::Dynamic> x = A.jacobiSvd(
        //Only compute info needed for least squares solution
        Eigen::ComputeThinU | Eigen::ComputeThinV)
        .solve(b);

    //Transfer the results into a CTF
    std::vector<CTF::ctf_t> results;
    results.resize(n);
    for(int i = 0; i < n; i++){
        results[i] = exp(x(i,0));
    }
    CTF ctf(results);

    //Potentially extract the pixel values to verify quality of fit
    const bool extractPoints = retPixels != NULL;
    if(extractPoints){
        assert(retPixels->empty());

        for(size_t i = n; i < n + numSamples; i++){
           PixelResult curr;
           curr.x = samplePositions[i-n].x;
           curr.y = samplePositions[i-n].y;
           curr.irradiance = exp(x(i,0));
           retPixels->push_back(curr);
        }
    }

    //All done
    return ctf;
}


std::ostream& operator<<(std::ostream& os, const CTFSolver& s){
    os << "CTFSolver{ lambda = " << s.lambda << ", channel = " <<
        s.chan << ", numSamples = " << s.numSamples << " }";
    return os;
}


