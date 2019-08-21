cmake_minimum_required(VERSION 3.12)
project(EASCPP)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp func.cpp func.h base.cpp base.h consts.h)
add_executable(EASCPP ${SOURCE_FILES})