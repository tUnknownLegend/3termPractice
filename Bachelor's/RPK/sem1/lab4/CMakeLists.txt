cmake_minimum_required (VERSION 3.2)

project ("LAB4")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}" )

include_directories(./src ./include)

file(GLOB SOURCES *.cpp *.h)

set(SOURCES ./main.cpp)
#set(HEADERS )

add_executable(${PROJECT_NAME} ${SOURCES})
