cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

project(zmath LANGUAGES CXX)

file(GLOB srcs CONFIGURE_DEPENDS test/*.cpp include/*.h)

add_executable(test ${srcs})

target_include_directories(test PUBLIC include)