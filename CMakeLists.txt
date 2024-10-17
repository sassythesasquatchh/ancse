cmake_minimum_required(VERSION 3.5)
project(cpphil)

# Find Eigen
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

if(NOT EIGEN3_INCLUDE_DIR)
  message(FATAL_ERROR "Eigen not found. Please install it using sudo apt-get install libeigen3-dev (Debian).")
endif()

# Find Python 2.7 and its components
find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)

if(NOT Python3_NumPy_INCLUDE_DIRS)
  message(FATAL_ERROR "NumPy headers not found for Python 2.7. Please install it using sudo apt-get install python2.7 python2.7-dev (Debian) or use an appropriate conda environment.")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(Examples)
add_subdirectory(Library)

