#Main CMAKE file
#Program to recover a camera transfer function.
#Currently only works on 8 bit cameras; this could however be generalized to N-bit sensors.
#
#Implements a slightly modified version of the 1997 paper
#"Recovering High Dynamic Range Radiance Maps from Photographs"
#by Debevec and Malik
#
#Colin Braley
#Initial revision                                        - May 15 2012
#Added code to make HDR images(replaces my old HDR code) - May 21 2012

#The basics
cmake_minimum_required(VERSION 2.6)
project(CTFSolve)
set(CMAKE_VERBOSE_MAKEFILE OFF)

#Application for finding camera CTF functions
set(CTF_SRCS  src/main.cpp src/CTFSolver.cpp src/CTF.cpp)
set(CTF_APP   bin/ctf_find )

#Application for making HDR images from exposure stacks
set(HDR_MAKE_SRCS src/mainHDRMake.cpp src/CTF.cpp src/CTFSolver.cpp src/WeightingFunctions.cpp)
set(HDR_MAKE_APP  bin/hdr_make )

#To build in debug mode change BUILD_TYPE to "Debug"
set(BUILD_TYPE  Release )

#Compiler flags
#TODO: The below code only works with G++
IF(CMAKE_COMPILER_IS_GNUCC)
    #set(CMAKE_CXX_FLAGS "-fpermissive") #Permissive only needed for Eigen library issues....

    #Always use -Wall and enable SSE 3 vector instructions
    set(CMAKE_CXX_FLAGS "-Wall -msse3")

    #In debug mode, don't optimize, output debugger symbols, and don't inline
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline")

    #In release mode, allow inlining, enable OpenMP and turn off assertions
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -finline-functions -fopenmp -D NDEBUG")
ELSE(CMAKE_COMPILER_IS_GNUCC)
    MESSAGE(FATAL_ERROR "This CMake file only works with g++!")
ENDIF(CMAKE_COMPILER_IS_GNUCC)

#Executables to create
set(CMAKE_BUILD_TYPE ${BUILD_TYPE})

#-------------------------------------------------------------
#Typically should not need to modify below this line----------
#-------------------------------------------------------------
add_executable(${CTF_APP}      ${CTF_SRCS}     )
add_executable(${HDR_MAKE_APP} ${HDR_MAKE_SRCS})

MESSAGE( STATUS "----------------------------------------")
MESSAGE( STATUS "\tBuild type: ${CMAKE_BUILD_TYPE}"       )
MESSAGE( STATUS "\tCXX Flags: ${CMAKE_CXX_FLAGS}"         )
MESSAGE( STATUS "----------------------------------------")




