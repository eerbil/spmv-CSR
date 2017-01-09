cmake_minimum_required(VERSION 2.8.8)
project(spMV-PlainCSR)

message(STATUS "System is " ${CMAKE_SYSTEM_NAME})
message(STATUS "Compiler is " ${CMAKE_CXX_COMPILER_ID})
message(${CMAKE_SOURCE_DIR})

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

message(STATUS "CXX Flags: " ${CMAKE_CXX_FLAGS})
message(STATUS "Linker Flags: " ${CMAKE_EXE_LINKER_FLAGS})

set(dir ${CMAKE_CURRENT_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${dir}")

set(SOURCE_FILES main.cpp
                 matrix.cpp)

add_executable(spMV ${SOURCE_FILES})

