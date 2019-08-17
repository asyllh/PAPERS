cmake_minimum_required(VERSION 3.12)
project(CMSASCPP)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp func.cpp func.h base.cpp base.h consts.h dlx.cpp dlx.h)
add_executable(CMSASCPP ${SOURCE_FILES})