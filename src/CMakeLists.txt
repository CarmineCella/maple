cmake_minimum_required(VERSION 2.8.9)
project (mptypes)

SET(CMAKE_CXX_FLAGS "-Wall -g -O2 -std=c++11")

option(ENABLE_READLINE "Enable readline shell" OFF)
if (ENABLE_READLINE) 
    set (LIBS ${LIBS} readline) 
    add_definitions (-DENABLE_READLINE)
endif()

add_executable(mptypes mptypes.cpp algorithms.h fourier.h WavFile.h Oscillator.h)
target_link_libraries (mptypes dl ${LIBS})
