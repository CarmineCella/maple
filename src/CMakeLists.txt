cmake_minimum_required(VERSION 2.8.9)
project (maple)

SET(CMAKE_CXX_FLAGS "-Wall -g2 -Og -std=c++11")

add_executable(maple maple.cpp utils.h algorithms.h fourier.h WavFile.h Matrix.h)
target_link_libraries (maple dl ${LIBS})
