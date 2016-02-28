cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project("DetectorSim" CXX)
set(DetectorSim_VERSION_MAJOR 1)
set(DetectorSim_VERSION_MINOR 0)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/library_cmake_files)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Boost COMPONENTS program_options filesystem regex system thread REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})
link_directories(${GSL_LIBRARY_DIRS})

add_executable(DetectorSim main.cpp Response.cpp Detector.cpp CalcCoordinator.cpp Extrap1d.cpp OutputResult.cpp)

set(needed_features cxx_strong_enums cxx_constexpr)
target_compile_features(DetectorSim PRIVATE ${needed_features})

target_link_libraries(DetectorSim ${Boost_LIBRARIES})
target_link_libraries(DetectorSim gsl)
