cmake_minimum_required (VERSION 3.2)

project ("CW_SEM_6")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -O3")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

include_directories(./src ./include)

file(GLOB SOURCES *.cpp *.h)

set(SOURCES ./src/main.cpp ./src/shared.cpp ./src/methods.cpp ./src/nonLinearSolve.cpp)
set(HEADERS ./include/SOLE.h ./include/shared.h ./include/methods.h ./include/nonLinearSolve.h ./include/SOLE.h)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
