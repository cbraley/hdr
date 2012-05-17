#Main CMAKE file
#Program to recover a camera transfer function.
#Currently only works on 8 bit cameras; this could however be generalized to N-bit sensors.
#
#Implements a slightly modified version of the 1997 paper 
#"Recovering High Dynamic Range Radiance Maps from Photographs"
#by Debevec and Malik
#
#Colin Braley
#May 15 2012

#The basics
cmake_minimum_required(VERSION 2.6)
project(CTFSolve)
set(CMAKE_VERBOSE_MAKEFILE OFF)

#Main variables
set(HDR_SRCS  src/main.cpp src/CTFSolver.cpp src/CTF.cpp)
set(HDR_APP   bin/ctf_find )

#To build in debug mode change BUILD_TYPE to "Debug"
set(BUILD_TYPE  Release )

#Compiler flags
#TODO: The below code only works with G++
IF(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS "-Wall")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-inline -msse3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -finline-functions -msse3")
ELSE(CMAKE_COMPILER_IS_GNUCC)
    MESSAGE(FATAL_ERROR "This CMake file only works with g++!")
ENDIF(CMAKE_COMPILER_IS_GNUCC)

#Executables to create
set(CMAKE_BUILD_TYPE ${BUILD_TYPE})

#-------------------------------------------------------------
#Typically should not need to modify below this line----------
#-------------------------------------------------------------

add_executable(${HDR_APP} ${HDR_SRCS})

MESSAGE( STATUS "----------------------------------------")
MESSAGE( STATUS "\tBuild type: ${CMAKE_BUILD_TYPE}"       )
MESSAGE( STATUS "\tCXX Flags: ${CMAKE_CXX_FLAGS}"         )
MESSAGE( STATUS "----------------------------------------")




