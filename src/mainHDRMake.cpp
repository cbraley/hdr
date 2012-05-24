#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
//--
#define cimg_display 0     //Don't compile cimg to use X11 display
#define cimg_verbosity 0   // Disable modal window in CImg exceptions.
#include "cimg/CImg.h"
#undef cimg_display
using namespace cimg_library;
//--
#include "Eigen/Eigen"
#include "Eigen/Dense"
//--
#include "CTF.h"
#include "CTFSolver.h"
//--
#include "LinearRegression.h"

static const std::string DIR_SEP("/");

//Simple struct for a pixel coordinate
typedef struct PixelCoord{
    int x,y;
    PixelCoord(int xp = -1, int yp = -1) : x(xp), y(yp) {}
}PixelCoord;


/**
 *  Make an HDR and return the # of bad pixels.
 *
 *  @param images is the exposure stack.
 *  @param pixelsToConsider is a list of pixels that should be considered.
 *   Pixels not in this list are left untouched in outHDR.
 *  @param ctf is the tabulated camera transfer function.
 *  @param outHDR is the output image.  This must be alloacted to proper size by
 *   the callee.
 *  @param outN is a pointer to an 8 bit LDR image to which we will output the number of valid
 *   measurments at each pixel.  If this is NULL, we won't consider it.
 *  @param outR is a pointer to an HDR image to which we will output the quality of fit per pixel.
 *   If this is NULL, we won't consider it.
 */
int makeHDR(const std::vector<CTFSolver::ImageExposurePair>& images,
    const std::vector<PixelCoord>& pixelsToConsider,
    const CTF& ctf,
    unsigned char validBegin, unsigned char validEnd,
    CImg<float>& outHDR,
    CImg<unsigned char>* outN = NULL,
    CImg<float>* outR = NULL
    )
{
    //Load images
    //guaranteed to succeed(barring that the user hasn't deleted a file a few ms ago)
    //since we verified their existance in main(...)
    std::vector< CImg<float> > ims;
    for(size_t j = 0; j < images.size(); j++){
        ims.push_back( CImg<float>(images[j].imagePath.c_str()) );
    }

    int badPixCount = 0; //Count # of pixels with no samples

    //Sample the weighting function into a LUT
    CTF::ctf_t lut[256];
    WeightingFunctions::makeLUTHat(&(lut[0]),  validBegin, validEnd);

    //Loop over all pixels that we want to make HDR values for
    for(size_t i = 0; i < pixelsToConsider.size(); i++){
        const int x = pixelsToConsider[i].x;
        const int y = pixelsToConsider[i].y;

        //See equation 6 of the Debvec and Malik paper
        CTF::ctf_t numerator   = static_cast<CTF::ctf_t>(0.0);
        CTF::ctf_t denominator = static_cast<CTF::ctf_t>(0.0);
        int P = 0;
        for(size_t j = 0; j < images.size(); j++){
            const unsigned char pixelValue = ims[j](x,y);
            const CTF::ctf_t weight = lut[pixelValue];
            const CTF::ctf_t logExposureTime = static_cast<CTF::ctf_t>(
                log(images[j].getTime()) );

            const CTF::ctf_t ctfValue = ctf(pixelValue);
            const float denTerm = weight;
            const float numTerm = weight * (ctfValue - logExposureTime);

            //assert((float)logExposureTime <= (float)ctfValue); //This indicates a bad CTF
            assert( numTerm >= static_cast<CTF::ctf_t>(0.0f) ); //Also indicates a bad CTF

            numerator   += numTerm;
            denominator += denTerm;

            //Add 1 to P iff weight is > 0
            P += weight > 0 ? 1 : 0;
        }

        //Output final HDR value
        if(P != 0){ //Good estimate
            const float radianceEstimate = exp(
                static_cast<float>(numerator / denominator)
                );
            outHDR(x,y) = radianceEstimate;
        }else{  //Bad pixel!
            outHDR(x,y) = 0.0f;
            ++badPixCount;
        }


    }
    return badPixCount;
}


/// Same as above but optimized for the case of a linear CTF
int makeHDRLinear(const std::vector<CTFSolver::ImageExposurePair>& images,
    const std::vector<PixelCoord>& pixelsToConsider,
    unsigned char validBegin, unsigned char validEnd,
    CImg<float>& outHDR,
    CImg<unsigned char>* outN = NULL,
    CImg<float>* outR = NULL
    )
{
    assert(images.size() >= 2);

    //Load images
    //guaranteed to succeed(barring that the user hasn't deleted a file a few ms ago)
    //since we verified their existance in main(...)
    std::vector< CImg<float> > ims;
    for(size_t j = 0; j < images.size(); j++){
        ims.push_back( CImg<float>(images[j].imagePath.c_str()) );
    }

    //Declare count to return
    int badPixCount = 0; //Count # of pixels with no samples

    //Preallocate memory to avoid malloc's in inner loop
    //Allocate the matrix to the maximum possible size (what would occur if
    //all samples are valid)
    Eigen::Matrix<float, Eigen::Dynamic, 2> points((int)ims.size(), 2);

    //Loop over all pixels that we want to make HDR values for
    for(size_t i = 0; i < pixelsToConsider.size(); i++){
        const int x = pixelsToConsider[i].x;
        const int y = pixelsToConsider[i].y;
       
        //Loop over exposures
        int matRow = 0; //Keep track of what row we are on in "points"
        //also is the # of valid samples at this pixel
        for(size_t j = 0; j < images.size(); j++){

            //Get pixel value
            const unsigned char pixelValue = ims[j](x,y);
            if(pixelValue > validBegin && pixelValue < validEnd){
                //Sample is valid!
                points(matRow  , 0) = images[j].getTime();
                points(matRow++, 1) = static_cast<float>(pixelValue);
            }
        }

        //We need at least two points for a resonable radiance estimate
        //(you can fit an infinite # of lines to one point)
        float hdrVal = 0.0f;
        float residual = -1.0f;
        if(matRow >= 2){ //Enough samples
            LinearRegression::Line<float> hdrLine = 
                LinearRegression::linearRegression<float>(
                matRow,     //# of points
                &points,    //pointer to points
                outR == NULL ? NULL : &residual); 
                //Return the residual ONLY if we are saving residual images
                //note that residual computation slows this whole thing down!

                //Slope is HDR estimate
                hdrVal = hdrLine.m;

                //Negative slope indicates issues!
                if(hdrVal < 0.0f){badPixCount++;}

        }else{ //Not enough samples
            ++badPixCount;
        }

        //Output to the HDR image
        outHDR(x,y) = hdrVal;

        //Potentially output to visualizations
        if(outN != NULL){
            assert(matRow < 256);
            (*outN)(x,y) = (unsigned char)matRow;
        }
        if(outR != NULL){
            (*outR)(x,y) = residual;
        }

    }

    return badPixCount;
}





void test(){

    //Enter N points
    std::cout << "Linear regresion - Enter N: " << std::endl;
    int N = 3;
    std::cin >> N;
    std::cout << std::endl;

    //Enter data points
    Eigen::Matrix<float,Eigen::Dynamic,2> data = 
        Eigen::Matrix<float,Eigen::Dynamic,2>(N,2);
//    Eigen::Vector2f*  data     = new Eigen::Vector2f[N];
    for(int i = 0; i < N; i++){
        float x,y; x = y = 0.0f;
        std::cout << "Enter x y: ";
        std::cin >> x >> y;
        std::cout << std::endl;

        data(i,0) = x;
        data(i,1) = y;
    }

    //Compute line
    float residual = 0.0f;
    LinearRegression::Line<float> line = LinearRegression::linearRegression<float>(
        N, &data, &residual);

    std::cout << "Line is: " << line << std::endl;
    std::cout << "residual is: " << residual << std::endl;





    CTF::ctf_t LUT[256];


    while(true){
        int a,b; a = b = 0;
        std::cin >> a >> b;
        std::cout << std::endl;
        WeightingFunctions::makeLUTHat(&(LUT[0]),  (unsigned char)a, (unsigned char)b);
        for(int i = 0; i < 256; i++){
            std::cout << i << "\t" << LUT[i] << std::endl;
        }
    }

    exit(0);
}


int main(int argc, char** argv){

    //Make std::vector of arguments
    std::string appName(argv[0]);
    std::vector<std::string> args;
    for(int i = 1; i < argc; i++){
        args.push_back(argv[i]);
    }

    //Check for improper arguments
    if( args.size() == 0 ||
        (args.size() < 7 && args[0] != "--help"))
    {
        std::cerr << "Error, invalid arguments!  See: " << appName << " --help for usage info." << std::endl;
        return 1;
    }

    //Check if we should print usage info
    if(args[0] == "--help"){
        std::cout << "Usage: " << std::endl <<
            "" << appName << " [OPTIONS] strategy in_folder_path out_file [FILE_LIST]" << std::endl;
        std::cout << "Required arguments: " << std::endl <<
            "\tstrategy can be either -ctf_linear or --ctf_tabular ctf_file" << std::endl <<
            "\t\t-ctf_linear assumes a linear camera transfer function." << std::endl <<
            "\t\t--ctf_tabular uses a non-linear CTF provided in file \"ctf_file\"." << std::endl <<
            "\tin_folder_path is a path to a folder of LDR images." << std::endl <<
            "\tout_file is a path(including extension) to a .pfm image" << std::endl <<
            "\t[FILE_LIST] is a list LDR images of the form path_1 exp_time_1 ... path_N exp_time_N." << std::endl <<
            "\t\tAll images in FILE_LIST should reside in \"in_folder_path\"" << std::endl <<
            "\t\tAt least 2 images must be present(N>=2)" << std::endl <<
            "\t\tExposure times are parsed as type \"long,\" so they should be integral." << std::endl;
        std::cout << "Optional arguments: " << std::endl <<
            "\t--matte path        - Use LDR image \"path\" as a matte.  Non-white pixels in the matte are ignored." << std::endl <<
            "\t--toe_size X        - Don't include pixel values in the range [0,X] in the fit." << std::endl <<
            "\t--shoulder_size X   - Don't include pixel values in the range [255-X,255] in the fit." << std::endl <<
            "\t--out_r path        - Write image of residual to file \"path\"." << std::endl <<
            "\t--out_n path        - Write image of number of valid pixels to \"path\"" << std::endl <<
            "\t-silent             - Don't print information to stdout." << std::endl <<
            "\t-discard_bloom_pix  - Discard pixels with an immediate neighbor that is in the range [250,255]" << std::endl <<
            "";
        return 0;
    }

    //Parse options
    std::string inFolderPath; //Path to input data
    std::string outFilePath;  //Path to output image
    //Pixels in the range [validPixEnd, validPixEnd] are used for HDR
    int validPixBegin = 0;
    int validPixEnd   = 255;
    int bloomStart    = 256; //Pixels are considered blooming if they are > bloomStart
    //Write residual image to outRPath, to write no image set to ""
    std::string outRPath = "";
    //Write N samples image to outNPath, to write no image set to ""
    std::string outNPath = "";
    std::string matteImagePath = ""; //Path to matte image, for no matte set to ""
    bool ctfLinear = false; //Should we assume a linear CTF?
    std::string ctfFile = ""; //tabulated CTF file if we are assuming non-linear CTF
    bool silent = false; //Should we keep quiet?
    unsigned char toeSize = 0;
    unsigned char shoulderSize = 0;
    CTF ctf; //Camera transfer function

    //First parse optional options
    size_t index = 0;
    while(index < args.size()){
        std::string arg = args[index++];
        if(arg       == "--toe_size"        ){
            int toeSizei = atoi(args[index++].c_str());
            if(toeSizei < 0 || toeSizei > 255){
                std::cerr << "Out of range toe size: " << toeSizei << std::endl;
                return 5;
            }
            toeSize = (unsigned char)toeSizei;
            validPixBegin += toeSize;
        }else if(arg == "--shoulder_size"   ){
            int shoulderSizei = atoi(args[index++].c_str());
            if(shoulderSizei < 0 || shoulderSizei > 255){
                std::cerr << "Out of range shoulder size: " << shoulderSizei << std::endl;
                return 5;
            }
            shoulderSize = (unsigned char)shoulderSizei;
            validPixEnd -= shoulderSize;
        }else if(arg == "--matte"           ){
            matteImagePath = args[index++];
        }else if(arg == "--out_n"           ){
            outNPath = args[index++];
        }else if(arg == "--out_r"           ){
            outRPath = args[index++];
        }else if(arg == "-discard_bloom_pix"){
            bloomStart -= 6; //This puts bloomStart at 250
        }else if(arg == "-silent"){
            silent = true;
        }else{ //Done parsing optional arguments
            --index;
            break;
        }
    }

    //Parse the required arguments
    const int argsLeft = args.size() - index;
    if(argsLeft < 7){
        std::cerr << "After parsing optional arguments, only " << argsLeft << " arguments remained." << std::endl;
        std::cerr << "This is an insufficient number of arguments." << std::endl;
        return 1;
    }
    std::string ctfStratString = args[index++];
    if(ctfStratString == "-ctf_linear"){
        ctfLinear = true;

    }else if(ctfStratString == "--ctf_tabular"){
        ctfLinear = false;
        ctfFile = args[index++];
    }else{
        std::cerr << "Invalid CTF strategy: \"" << ctfStratString << "\"" << std::endl;
        return 2;
    }
    inFolderPath = args[index++];
    outFilePath  = args[index++];

    //Parse the image exposure pairs
    std::vector<CTFSolver::ImageExposurePair> images;
    while(index < args.size()){
        //Make sure we have a full pair to parse
        if(index + 2 > args.size()){
            std::cerr << "Invalid image exposure pair: " << args[index] << std::endl;
            return 3;
        }
        std::string path  = inFolderPath + DIR_SEP + args[index++];
        long microseconds = atol(args[index++].c_str());
        images.push_back(CTFSolver::ImageExposurePair(microseconds, path));
    }


    //Make sure we have at least two images
    if(images.size() < 2){
        std::cerr << "Error - At least 2 images are required!" << std::endl;
        return 3;
    }
    //Sort the images by exposure time
    std::sort(images.begin(), images.end());

    //Print info as long as we are not in silent mode
    if(!silent){
        std::cout << "Command line option summary: " << std::endl;
        std::cout << "\tInput folder: " << inFolderPath << std::endl;
        std::cout << "\tOutput HDR: " << outFilePath << std::endl;
        std::cout << "\tValid pixel range [" << validPixBegin << ", " << validPixEnd << "]" << std::endl;
        if(bloomStart < 256){
            std::cout << "\tIgnoring bloom pixels and neighbors(bloom is > " << bloomStart << ")" << std::endl;
        }else{
            std::cout << "\tNot compensating for bloom." << std::endl;
        }
        std::cout << "\tCTF is: ";
        if(ctfLinear){
            std::cout << "assumed to be linear." << std::endl;
        }else{
            std::cout << ctfFile << std::endl;
        }
        if(outRPath != ""){
            std::cout << "\tOutputting residual to: " << outRPath << std::endl;
        }
        if(outNPath != ""){
            std::cout << "\tOutputting num samples visualization to: " << outNPath << std::endl;
        }
        if(matteImagePath != ""){
            std::cout << "\tMatte image: " << matteImagePath << std::endl;
        }else{
            std::cout << "\tNot using a matte image." << std::endl;
        }
        std::cout << "\t" << images.size() << " images: ";
        for(size_t i = 0; i < images.size(); i++){
            std::cout << images[i] << " ";
        }
        std::cout << std::endl;
    }


    //Make sure all the images load
    //also, get the width and height from disk
    int width, height, numChans; width = height = numChans = -1;
    std::string errStr;
    const bool imagesOK = CTFSolver::checkImagesOK(images,
        width, height, numChans, &errStr);
    if(!imagesOK){
        std::cerr << "Could not load 1 or more images!" << std::endl;
        std::cerr << "The issue was: \"" << errStr << "\"" << std::endl;
        return 7;
    }else if(numChans != 1){
        std::cerr << "Error - Only works on monochrome images!" << std::endl;
        return 4;
    }

    //Make the CTF
    if(ctfLinear){

        //Compute max( natural_log(exposure_times) )
        //(see section 2.2 of Debevec paper)
        const CTF::ctf_t maxExpTime = images[images.size()-1].getTime();
        ctf = CTF::makeLinearCTF( log(maxExpTime) );

    }else{
        //Load the CTF from a file
        const bool loadCTFOK = CTF::loadCTF(ctf, ctfFile);
        if(!loadCTFOK){
            std::cerr << "Could not load CTF from: " << ctfFile << std::endl;
            return 33;
        }

    }


    //If we got here we are able to load all the image
    //Lets make an HDR

    //Declare mem for output image
    CImg<float> hdr(width, height, 1, 1);

    //Zero out all pixels
    hdr.fill(0.0f);

    //Make a list of the pixels we want to compute HDRs at
    std::vector<PixelCoord> pixelsToConsider;
    if(matteImagePath != ""){
        try{
            CImg<unsigned char> matte(matteImagePath.c_str());
            if(matte.width() != width || matte.height() != height || matte.spectrum() != 3){
                std::cerr << "Invalid matte dimensions!" << std::endl;
                return 9;
            }
            pixelsToConsider.reserve(width * height);
            for(int x = 0; x < width; x++){
                for(int y = 0; y < height; y++){
                    //White pixels are included
                    if( matte(x,y,0,0) == 255 &&
                        matte(x,y,0,1) == 255 &&
                        matte(x,y,0,2) == 255 )
                    {
                        pixelsToConsider.push_back( PixelCoord(x,y) );
                    }
                }
            }
        }catch(const CImgException& ex){
            std::cerr << "Could not load matte image: " << matteImagePath << std::endl;
            return 10;
        }catch(...){ //Catches ANYTHING thown
            std::cerr << "Could not load matte image: " << matteImagePath << std::endl;
            return 10;
        }
    }else{
        //Consider all pixels!
        pixelsToConsider.reserve(width * height);
        for(int x = 0; x < width; x++){
            for(int y = 0; y < height; y++){
                pixelsToConsider.push_back( PixelCoord(x,y) );
            }
        }
    }

    //Make sure at least 1 pixel is on
    if(pixelsToConsider.size() < 1){
        std::cerr << "Error - No pixels were on in the matte!" << std::endl;
        return 11;
    }



    //Make an HDR

    //Declare output images for the "N" visualization as well as the "r" visualization
    //However, the pointers passed to the makeHDR function are non-null if and only if
    //the paths are not the empty string
    CImg<unsigned char> outN(width,height,1,1);
    outN.fill(0);
    CImg<float> outR(width,height,1,1);
    outR.fill(0.0f);
    CImg<unsigned char>* outNPtr = outNPath == "" ? NULL : &outN;
    CImg<float>* outRPtr = outRPath == "" ? NULL : &outR;
    int numCompleteErrors = -1;
    assert(validPixBegin < validPixEnd);
    if(ctfLinear){ //Linear CTF special case (faster)
        numCompleteErrors = makeHDRLinear(images, pixelsToConsider,
            validPixBegin, validPixEnd,
            hdr,
            outNPtr, outRPtr);
    }else{ //Non-linear CTF general case
        numCompleteErrors = makeHDR(images, pixelsToConsider,
            ctf,
            validPixBegin, validPixEnd,
            hdr,
            outNPtr, outRPtr);
    }
    if(numCompleteErrors > 0){
        std::cerr << "ERROR - Found: " << numCompleteErrors <<
            " error pixels when making HDR(s)!" << std::endl;
        std::cerr << "\tThis means that: " << numCompleteErrors << " pixel locations had < 2 images with pixels " << 
            " in range [" << validPixBegin << ", " << validPixEnd << "]" << std::endl;
        const float perc = (((float)numCompleteErrors)/((float)pixelsToConsider.size()) ) * 100.0f;
        std::cerr << "\t" << perc << " percent of the pixels are therefore invalid!" << std::endl;
    }




    //Write the output image(s)
    try{
        hdr.save(outFilePath.c_str());
        if(!silent){
            std::cout << "Wrote HDR result to: " << outFilePath << std::endl;
        }
        if(outNPtr != NULL){
            outNPtr->save(outNPath.c_str());
            if(!silent){
                std::cout << "Wrote N-samples visualization to: " << outNPath << std::endl;
            }
        }
        if(outRPtr != NULL){
            outRPtr->save(outRPath.c_str());
            if(!silent){
                std::cout << "Wrote residual visualization to: " << outRPath << std::endl;
            }
        }
    }catch(const CImgException& ex){
        std::cerr << "Could not save 1 or more of the output images!" << std::endl;
        return 17;
    }catch(...){ //Catches ANYTHING thown
        std::cerr << "Could not save 1 or more of the output images!" << std::endl;
        return 18;
    }

    return 0;
}

