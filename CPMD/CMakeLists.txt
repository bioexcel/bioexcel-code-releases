cmake_minimum_required(VERSION 2.6)
project(MimicCommLib)

option(INCLUDE_TESTS "Choose whether to build tests in this run or not" OFF)

set(CMAKE_CXX_STANDARD 11)

find_package(MPI REQUIRED)
include_directories("${MPI_INCLUDE_PATH}")

add_subdirectory(src)

if(INCLUDE_TESTS)
    enable_testing()
    set(CMAKE_CXX_FLAGS "-std=c++11")

    find_package(Boost 1.50.0  REQUIRED COMPONENTS
            serialization system)

    if(NOT Boost_FOUND)
        message(FATAL_ERROR "Could not find boost!")
    endif()
    include_directories(${Boost_INCLUDE_DIRS})
    set(TEST_SOURCE_FILES ${SOURCE_FILES} tests/tests.cpp)
    add_executable(test_exec ${TEST_SOURCE_FILES})
    set (test_parameters -n 1 ${CMAKE_CURRENT_BINARY_DIR}/test_exec)
    target_link_libraries(test_exec ${Boost_LIBRARIES})
    add_test(NAME "MyTest" COMMAND "mpirun" ${test_parameters})
endif(INCLUDE_TESTS)