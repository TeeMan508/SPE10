cmake_minimum_required(VERSION 2.14)

project(SPE10 DESCRIPTION "Solving differential equation" LANGUAGES CXX)
include_directories(src)
add_executable(
    main main.cpp 
)