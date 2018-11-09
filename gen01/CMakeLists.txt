cmake_minimum_required (VERSION 3.10)

project (gen01)

include_directories(include)

file(GLOB SOURCES "src/*.cpp")

add_executable(gen01 ${SOURCES})

set(CMAKE_BUILD_TYPE Debug)
