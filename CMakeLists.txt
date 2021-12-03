# Set the minimum required version of cmake for a project.
cmake_minimum_required(VERSION 3.10)

# Set Variables
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 99)
set(PROJECT_NAME "week2")
set(LIBIGL_INCLUDE_DIR /home/sungho/data/libigl/v2.2.0/include)

# Project
PROJECT(${PROJECT_NAME} CXX)
# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
find_package(LIBIGL REQUIRED QUIET)
include_directories(${PROJECT_SOURCE_DIR}/include
	   	    ${LIBIGL_INCLUDE_DIR}
                    ${LIBIGL_INCLUDE_DIR}/../external
                    ${LIBIGL_INCLUDE_DIR}/../external/glad/include
		    )
#link_directories(/home/sungho/data/libigl/v2.2.0/lib/)
		    
# The main program
add_executable(${PROJECT_NAME} main.cpp)

# Link
target_link_libraries( ${PROJECT_NAME} 
	igl::core
	igl::opengl
	igl::opengl_glfw
)
