cmake_minimum_required(VERSION 3.3)
project(dm818_serial)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES common.cpp common.h grid.cpp grid.h serial_base.cpp serial_base.h)

add_executable(dm818_base_serial ${SOURCE_FILES} serial_base_main.cpp)
add_executable(dm818_serial serial.cpp serial.h serial_main.cpp ${SOURCE_FILES} )
add_executable(dm818_test ${SOURCE_FILES} test_suite.cpp)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_executable(mpi mpi.cpp mpi.h ${SOURCE_FILES})
target_link_libraries(mpi ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(mpi PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(mpi PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

target_compile_definitions(mpi PRIVATE DEBUG)
target_compile_definitions(dm818_serial PRIVATE DEBUG)
