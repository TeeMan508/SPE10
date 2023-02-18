cmake_minimum_required(VERSION 2.14)

project(SPE10 DESCRIPTION "Solving differential equation" LANGUAGES CXX)


add_executable(spe10 source/main.cpp)
add_executable(test testing/tests.cpp)

include_directories(${PROJECT_SOURCE_DIR})