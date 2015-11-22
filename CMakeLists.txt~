cmake_minimum_required(VERSION 2.8)
project(dm818_serial)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES_MPI mpi-other/common.cpp mpi-other/common.h mpi-other/grid.cpp mpi-other/grid.h)
set(SOURCE_FILES_SERIAL serial-other/common.cpp serial-other/common.h serial-other/grid.cpp serial-other/grid.h serial-other/serial_base.cpp serial-other/serial_base.h)

add_executable(dm818_serial serial-other/serial.cpp serial-other/serial.h serial-other/serial_main.cpp ${SOURCE_FILES_SERIAL} )
add_executable(dm818_base_serial serial-other/serial_base_main.cpp ${SOURCE_FILES_SERIAL} )
add_executable(dm818_test ${SOURCE_FILES_SERIAL} serial-other/test_suite.cpp)


find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_executable(mpi mpi-other/mpi.cpp mpi-other/mpi.h ${SOURCE_FILES_MPI})
target_link_libraries(mpi ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(mpi PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(mpi PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

#target_compile_definitions(mpi PRIVATE DEBUG)
#target_compile_definitions(dm818_serial PRIVATE DEBUG)
