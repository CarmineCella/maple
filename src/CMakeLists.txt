cmake_minimum_required(VERSION 2.8.9)
project (mptypes)

SET(CMAKE_CXX_FLAGS "-Wall -g2 -O2 -std=c++11")

option(ENABLE_READLINE "Enable readline shell" OFF)
if (ENABLE_READLINE) 
    set (LIBS ${LIBS} readline) 
    add_definitions (-DENABLE_READLINE)
endif()

add_executable(mptypes mptypes.cpp utils.h algorithms.h fourier.h WavFile.h Matrix.h)
target_link_libraries (mptypes dl ${LIBS})
