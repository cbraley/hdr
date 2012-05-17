#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstring>
#include <fstream>
//--
#include "CTF.h"
#include "CTFSolver.h"

static const int DFLT_NUM_SAMPS = 500;
static const int DFLT_CHAN      = 0  ;
static const double DFLT_LAMBDA = 3.0;

//app [OPTIONS] --num_files N  file_1 time_1 .... file_N time_N
//
//Valid options
//    --help
//    --num_samps INT
//    --lambda FLOAT
//    --weighting_func  {hat,uniform}
//    --silent
int main(int argc, char** argv){

    //Check for help message
    if(argc > 1 && strcmp(argv[1],"--help") ==  0){
        std::cout << "Usage: " << std::endl << 
            argv[0] << " [OPTIONS] --num_files N  file_1 time_1 ... file_N time_N" << std::endl;
        std::cout << "OPTIONS include: "        << std::endl;
        std::cout << "\t--help"                 << std::endl;
        std::cout << "\t\tPrint this usage information." << std::endl;
        std::cout << "\t--num_samps INTEGER"    << std::endl;
        std::cout << "\t\tNumber of image samples to take.  Defaults to " << DFLT_NUM_SAMPS << std::endl;
        std::cout << "\t--bloom_discard INTEGER" << std::endl;
        std::cout << "\t--lambda FLOAT"         << std::endl;
        std::cout << "\t\tSmoothing coeffecient.  Defaults to " << DFLT_LAMBDA << std::endl;
        std::cout << "\t--weight_func {hat, hat_10}" << std::endl;
        std::cout << "\t\tDefaults to \"hat,\" the function used in Debevec and Malik 1997." << std::endl;
        std::cout << "\t\that uses a triangle filter that starts at 0 and ends at 255." << std::endl;
        std::cout << "\t\that_10 w uses with 0 weight on the upper and lower 10 values." << std::endl;
        std::cout << "\t--out_file fileName"    << std::endl;
        std::cout << "\t\tFile to write CTF data to.  Defaults to writing to stdout if this is not specified" << std::endl;
        std::cout << "\t--out_file_points fileName"    << std::endl;
        std::cout << "\t\tFile to write raw points used in solve to." << std::endl;
        std::cout << "\t--silent" << std::endl;
        std::cout << "\t\tIf specified, we only write(or print) the CTF and do nothing else." << std::endl;
        return 0;
    }

    //Check for improper arg count
    if(argc < 5){
        std::cerr << "Invalid number of arguments." << std::endl;
        std::cerr << "See " << argv[0] << " --help" << std::endl;
        return 1;
    }

    //First parse optional arguments
    int index = 1;
    int numSamps = DFLT_NUM_SAMPS;
    double lambda = DFLT_LAMBDA;
    int chan = DFLT_CHAN;
    bool silent = false;
    std::string outFile("-");
    std::string outFilePoints("");
    CTFSolver::WeightingFunc wFunc = CTFSolver::HAT;
    while(strcmp(argv[index],"--num_files") != 0 && index < argc){
        char* arg = argv[index++];
        if(strcmp(arg,"--num_samps")    == 0){
            numSamps = atoi(argv[index++]);
        }else if(strcmp(arg,"--lambda") == 0){
            lambda = strtod(argv[index++], NULL);
        }else if(strcmp(arg,"--weight_func") == 0){
            char* wFuncName = argv[index++];
            if(strcmp(wFuncName, "hat") == 0){
                wFunc = CTFSolver::HAT;
            }else if(strcmp(wFuncName,"hat_10") == 0){
                wFunc = CTFSolver::HAT_10;
            }else{
                std::cerr << "Unknown weighting function: " << wFuncName << std::endl;
                return 2;
            }
        }else if(strcmp(arg,"--silent") == 0){
            silent = true;
        }else if(strcmp(arg,"--out_file") == 0){
            char* f = argv[index++];
            outFile = std::string(f);
        }else if(strcmp(arg,"--out_file_points") == 0){
            char* f = argv[index++];
            outFilePoints = std::string(f);
        }else{
            std::cerr << "Unknown option: " << arg << std::endl;
            return 1;
        }
    }
    const bool writeCurveToStdOut = outFile == "-";
    const bool writePointsToFile = outFilePoints != "";

    //Parse number of files
    int numFiles = -1;
    //Make sure we just landed on the --num_files argument
    if(index >= argc || strcmp(argv[index],"--num_files") != 0){
        std::cerr << "Error - Missing --num_files argument." << std::endl;
        std::cerr << "See " << std::endl << 
            argv[0] << " --help" << std::endl;
        return 1;
    }
    numFiles = atoi(argv[++index]);
    if(numFiles < 2){
        std::cerr << "Error - At least 2 images required!" << std::endl;
        return 2;
    }
    ++index;


    //Read the images
    std::vector<CTFSolver::ImageExposurePair> images;
    int stopIndex = index + (numFiles*2);
    int filesRead = 0;
    while(index < stopIndex){

        //make sure we don't read too far
        if(index >= argc){
            std::cerr << "Error - Could not read all: " << numFiles <<
                " files from the command line." << std::endl;
            return 3;
        }

        //Read the pair
        std::string path(argv[index++]);
        double time = strtod(argv[index++], NULL);
        images.push_back(CTFSolver::ImageExposurePair(time,path));

        ++filesRead;
    }

    //Potentially print info
    if(!silent){
        std::cout << "Starting linear solve for CTF creation.  Parameters: " << std::endl;
        std::cout << "\tlambda      = " << lambda   << std::endl;
        std::cout << "\tnum_samples = " << numSamps << std::endl;
        std::cout << "\tchannel     = " << chan     << std::endl;
        if(writeCurveToStdOut){
            std::cout << "\tWriting curve to stdout." << std::endl;
        }else{
            std::cout << "\tWriting curve to " << outFile << std::endl;
        }
    }

    //Set up the solver
    CTFSolver solver(images, numSamps, lambda, chan);
    solver.setWeightingFunc(wFunc);

    //Setup parameters in the event that we want to return pixel irradiances
    std::vector<CTFSolver::PixelResult> retPixels;
    std::vector<CTFSolver::PixelResult>* retList = 
        writePointsToFile ? &retPixels : NULL;

    //Do the solve(takes some time)
    CTF ctf = solver.solve(retList); 

    //Output as desired
    if(writeCurveToStdOut){
        std::cout << ctf << std::endl;
    }else{
        std::fstream file(outFile.c_str(), std::fstream::out);
        if(file.good()){
            file << ctf << std::endl;
            file.close();
        }else{
            file.close();
            std::cerr << "Could not write to file: " << outFile << std::endl;
            return 3;
        }
    }

    //Potentially write points out
    if(writePointsToFile){
        std::fstream file(outFilePoints.c_str(), std::fstream::out);
        if(file.good()){
            const bool wok = solver.writePixelPoints(retPixels, file);
            file.close();
            if(!wok){
                std::cerr << "Could not write pixel points to file: " << outFilePoints << std::endl;
                return 4;
            }
        }else{
            file.close();
            std::cerr << "Could not write pixel points to file: " << outFilePoints << std::endl;
            return 5;
        }


    }

    //All done
    return 0;
}




