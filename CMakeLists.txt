# Set the minimum required version of cmake for a project.
cmake_minimum_required(VERSION 3.1)

# Set Variables
SET(CMAKE_CXX_FLAGS "-std=c++0x") # C++11
SET(PROJECT_NAME "eigen")

# Project
PROJECT(${PROJECT_NAME} CXX)

# Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)


# Include directories
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/include)

include_directories(${PROJECT_SOURCE_DIR}/dqs)
file(GLOB dqs ${PROJECT_SOURCE_DIR}/dqs/*)

# The main program
add_executable(${PROJECT_NAME} main.cpp ${SRC_FILES} ${HEADER_FILES} ${dqs})

# Link
target_link_libraries( ${PROJECT_NAME} )
