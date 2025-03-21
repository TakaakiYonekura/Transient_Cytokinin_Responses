cmake_minimum_required(VERSION 3.10)
project(FindSolution_Cytokinin)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(NLOPT_LIB_PATH "${CMAKE_SOURCE_DIR}/third_party/nlopt/build/Release")

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/third_party/eigen
    ${CMAKE_SOURCE_DIR}/third_party/nlopt/src/api
    ${CMAKE_SOURCE_DIR}/third_party/nlopt/build
)

# Create executable target
add_executable(FindSolution_Cytokinin
    src/main.cpp
    src/difAve.cpp
    src/DifferentialEquation.cpp
    src/readCSV.cpp
    src/stdafx.cpp  # PCH source
)

# for Windows
if (MSVC)
    message(STATUS "MSVC compiler detected: Applying precompiled header settings.")
    set(PCH_PATH "${CMAKE_BINARY_DIR}/stdafx.pch")
    message(STATUS "PCH Output Path: ${PCH_PATH}")
    set_source_files_properties(src/stdafx.cpp PROPERTIES COMPILE_FLAGS "/Ycstdafx.h /Fp${PCH_PATH}")
endif()

# Link NLopt library
target_link_libraries(FindSolution_Cytokinin ${NLOPT_LIBRARY})
target_link_libraries(FindSolution_Cytokinin PRIVATE "${NLOPT_LIB_PATH}/nlopt.lib")

add_custom_command(TARGET FindSolution_Cytokinin POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/third_party/nlopt/build/Release/nlopt.dll
    ${CMAKE_BINARY_DIR}/nlopt.dll
)
