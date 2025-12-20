cmake_minimum_required(VERSION 3.20)
project(mpiL1 C)

set(CMAKE_C_STANDARD 11)

find_package(MPI REQUIRED)

enable_testing()

add_executable(mpiL1
        L1.c
)

target_include_directories(mpiL1 PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(mpiL1 PRIVATE flint mpfr gmp m)

target_link_libraries(mpiL1 PRIVATE MPI::MPI_C)


