cmake_minimum_required(VERSION 3.23)
project(manipulatorTools)

find_package(Eigen3 REQUIRED)

set(CMAKE_CXX_STANDARD 11)

include_directories(
        ${EIGEN3_INCLUDE_DIRS}
)

add_executable(manipulatorTools main.cpp)
