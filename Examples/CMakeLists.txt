cmake_minimum_required(VERSION 3.5)

#project(CppHilbertProject LANGUAGES CXX)

# Set the C++ standard for all targets
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define a list of executable names and their source files
set(EXECUTABLES
    cmtest:test_build.cpp
    capacitor_Dirichlet:capacitor_Dirichlet.cpp
    capacitor_Neumann:capacitor_Neumann.cpp
    manufactured_Dirichlet:manufactured_Dirichlet.cpp
    capacitor_pointerror:capacitor_pointerror.cpp
    capacitor_indirect:capacitor_indirect.cpp
    mesh_capacitor:mesh_capacitor.cpp
    Lshape_Dirichlet:Lshape_Dirichlet.cpp
    Lshape_Neumann:Lshape_Neumann.cpp
    Lshape_pointerror:Lshape_pointerror.cpp
    lshape_adaptive_plot:lshape_adaptive_plot.cpp
    meshRefineTest:meshRefineTest.cpp
    Lshape_PCG:Lshape_PCG.cpp
)

# Function to add executables and set properties
function(add_hilbert_executable target source)
    add_executable(${target} ${source})
    target_include_directories(${target} PRIVATE ${EIGEN3_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/Library/source ${PROJECT_SOURCE_DIR}/Library/ ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
    target_link_libraries(${target} PRIVATE CppHilbert CppHilbertMySolution ${Python3_LIBRARIES})
endfunction()

# Loop over the executables and create them
foreach(executable_pair ${EXECUTABLES})
    string(REPLACE ":" ";" executable_item ${executable_pair})
    list(GET executable_item 0 executable_name)
    list(GET executable_item 1 source_file)
    add_hilbert_executable(${executable_name} ${source_file})
endforeach()
