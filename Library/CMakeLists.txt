cmake_minimum_required(VERSION 3.5)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${PROJECT_SOURCE_DIR}/Library)
include_directories(${EIGEN3_INCLUDE_DIR})

file(GLOB SOURCES "./*.hpp" "source/*.cpp" "source/*.hpp" "source/gaussTable.h" "source/gaussQuadrature.h" "source/gaussQuadrature.c")
file(GLOB SOURCESMYSOLUTION "templates/*.cpp")

add_library(CppHilbert SHARED ${SOURCES})
add_library(CppHilbertMySolution SHARED ${SOURCESMYSOLUTION})

# Link CppHilbert to CppHilbertMySolution
target_link_libraries(CppHilbertMySolution PUBLIC CppHilbert)

# Define installation location
if(NOT CPPHILBERT_ROOT AND NOT $ENV{CPPHILBERT_ROOT} STREQUAL "")
    set( CPPHILBERT_LIB   $ENV{CPPHILBERT_ROOT}/include  )
else()
    set( CPPHILBERT_LIB   ${CMAKE_BINARY_DIR}/../include )
endif()
install(TARGETS CppHilbert DESTINATION ${CPPHILBERT_LIB})

