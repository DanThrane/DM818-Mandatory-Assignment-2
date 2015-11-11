cmake_minimum_required(VERSION 3.3)
project(dm818_serial)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES common.cpp common.h grid.cpp grid.h serial.cpp serial.h serial_base.cpp serial_base.h)

add_executable(dm818_base_serial ${SOURCE_FILES} serial_base_main.cpp)
add_executable(dm818_serial ${SOURCE_FILES} serial_main.cpp)
add_executable(dm818_serial_debug ${SOURCE_FILES} serial_main.cpp)
add_executable(dm818_test ${SOURCE_FILES} test_suite.cpp)

target_compile_definitions(dm818_serial_debug PRIVATE DEBUG)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_executable(mpi mpi.cpp mpi.h ${SOURCE_FILES})
target_link_libraries(mpi ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(hello PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(hello PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()