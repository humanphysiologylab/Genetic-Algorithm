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


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -fopenmp -g -fPIC")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG -fopt-info-missed=missed.all") #-fopt-info")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wconversion")

list( INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake )

find_package(MPI REQUIRED)
if (MCMC_ENABLED)
  find_package(ROOT 6.16 CONFIG REQUIRED)
  find_package(BAT REQUIRED)
  include_directories(${BAT_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})
  message(STATUS "ROOT:" ${ROOT_INCLUDE_DIRS})
endif()

include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories("Eigen")
include_directories("pcg-cpp-0.98")
include_directories("LSODA")
include_directories("json")

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download"
)
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download"
)

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This adds the following targets:
# gtest, gtest_main, gmock and gmock_main
add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
                 "${CMAKE_BINARY_DIR}/googletest-build"
)

# The gtest/gmock targets carry header search path dependencies
# automatically when using CMake 2.8.11 or later. Otherwise we
# have to add them here ourselves.
if(CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include"
                        "${gmock_SOURCE_DIR}/include"
    )
endif()



set(SOURCE_FILES
        src/tournament_selection.cpp
        src/sbx_crossover.cpp
        src/cauchy_mutation.cpp
        src/writing_to_output_files.cpp
        src/maleckar.cpp
        src/optimization_problem.cpp
        src/gradient_descent.cpp
        src/polynomial_mutation.cpp
        src/basic_function_functor.cpp
        src/test_functions.cpp
        src/maleckar_model.cpp
        src/kernik_clancy_model.cpp
        src/mcmc.cpp
        LSODA/LSODA.cpp
        )

set(libname ${CMAKE_PROJECT_NAME}_lib)
add_library(${libname} STATIC ${SOURCE_FILES})
target_link_libraries(${libname} ${MPI_LIBRARIES} ${ROOT_LIBS} ${BAT_LIBRARIES})


add_executable(ga src/main.cpp)
target_link_libraries(ga ${libname})


enable_testing()
file(GLOB_RECURSE test_cases tests/*.cpp)
foreach(case_file ${test_cases})
    get_filename_component( case_name ${case_file} NAME_WE )
    #set (case_name test_${case_name})
    add_executable(${case_name} ${case_file})
    target_include_directories(${case_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(${case_name}
        ${libname} gtest gmock_main ${ROOT_LIBS} ${BAT_LIBRARIES})
    add_test(NAME ${case_name}
             COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${case_name}
             WORKING_DIRECTORY
             ${PROJECT_BINARY_DIR})
endforeach()
