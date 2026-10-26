cmake_minimum_required(VERSION 3.10)
project(nbfb CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_package(Threads REQUIRED)

add_library(nbfb_lib STATIC nbfb.cpp nbfb_simple.cpp)
target_include_directories(nbfb_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(nbfb_lib PRIVATE -Wall -Wextra)

add_executable(nbfb nbfb_cli.cpp)
target_link_libraries(nbfb PRIVATE nbfb_lib Threads::Threads)
target_compile_options(nbfb PRIVATE -Wall -Wextra)
