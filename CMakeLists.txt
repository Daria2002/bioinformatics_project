cmake_minimum_required(VERSION 3.5)
project(bioinformatics_project)
#define the compiler
set(CMAKE_CXX_FLAGS "-std=c++11")
 
# Add libbf as external project
include(ExternalProject)
ExternalProject_Add(project_bf
    GIT_REPOSITORY          "https://github.com/mavam/libbf.git"
    PREFIX                  "${CMAKE_BINARY_DIR}/libbf"
    SOURCE_DIR              "${CMAKE_BINARY_DIR}/libbf/src/project_bf"
    BINARY_DIR              "${CMAKE_BINARY_DIR}/libbf/src/project_bf-build"
    BUILD_COMMAND           "make"
    INSTALL_COMMAND         ""
    )
   
# TODO: Get these variables using ExternalProject_Get_Property macro
set(BF_INCLUDE_DIR "${CMAKE_BINARY_DIR}/libbf/src/project_bf")
set(BF_LIB_DIR "${CMAKE_BINARY_DIR}/libbf/src/project_bf-build/lib")
 
include_directories(
    include
    "${BF_INCLUDE_DIR}"
    )
link_directories(
    "${BF_LIB_DIR}"
    )
 
add_executable(main src/main.cpp)
add_dependencies(main project_bf)
target_link_libraries(main bf fasta)
add_library(fasta src/FastaParser.cpp)


