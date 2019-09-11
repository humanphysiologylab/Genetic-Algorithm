cmake_minimum_required(VERSION 2.8)

project(ga)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -O2")

find_package(MPI)

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
src/ord_model/atrium.cpp
)
add_executable(ga ${SOURCE_FILES})

target_link_libraries(ga ${MPI_LIBRARIES})
