cmake_minimum_required (VERSION 3.2)

project ("LAB1")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}" )

include_directories(./src ./include)

file(GLOB SOURCES *.cpp *.h)

set(SOURCES ./src/main.cpp ./src/shared.cpp ./src/Gauss.cpp src/QR.cpp)
set(HEADERS ./include/shared.h ./include/Gauss.h ./include/QR.h)

#find_package(GTest REQUIRED)
#include_directories(${GTEST_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

#install(TARGETS ${PROJECT_NAME} DESTINATION "LAB1/bin")