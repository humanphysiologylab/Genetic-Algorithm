cmake_minimum_required(VERSION 2.8)

project(ga)


set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

find_package(MPI REQUIRED)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(SOURCE_FILES
        src/main.cpp
        src/initial_population.cpp
        src/fitness_function.cpp
        src/tournament_selection.cpp
        src/sbx_crossover.cpp
        src/cauchy_mutation.cpp
        src/writing_to_output_files.cpp
        src/random_number_generator.cpp
        src/maleckar.cpp)

add_executable(ga ${SOURCE_FILES})

target_link_libraries(ga ${MPI_LIBRARIES})
