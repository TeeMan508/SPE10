cmake_minimum_required(VERSION 2.14)

project(SPE10 DESCRIPTION "Solving differential equation" LANGUAGES CXX)
include_directories(source)

add_executable(
    spe10 source/main.cpp 
)