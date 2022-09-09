cmake_minimum_required(VERSION 3.15)
project(TSP)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -O2")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}" )

file(GLOB SOURCES *.cpp *.h)
#find_package(GTest REQUIRED)
#include_directories(${GTEST_INCLUDE_DIRS})

set(SOURCES shared.cpp)

set(HEADERS shared.h)

add_executable(${PROJECT_NAME}_rec recursion.cpp ${SOURCES} ${HEADERS})
add_executable(${PROJECT_NAME}_simAnnealing simAnnealing.cpp ${SOURCES} ${HEADERS})

#add_executable(Tests ${SOURCES} ${HEADERS} tests.cpp)

#target_link_libraries(Tests GTest::Main)