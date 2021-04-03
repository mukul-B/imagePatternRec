cmake_minimum_required(VERSION 3.17)
project(rgb)

set(CMAKE_CXX_STANDARD 14)

add_executable(rgb image.h image.cpp rgb.h ReadImage.cpp ReadImageHeader.cpp WriteImage.cpp Threshold.cpp rgb.cpp MatricOP.h DiscriminantCases.h)