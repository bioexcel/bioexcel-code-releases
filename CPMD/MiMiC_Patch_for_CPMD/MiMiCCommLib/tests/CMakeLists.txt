find_package (Threads)
#include_directories(${GTEST_INCLUDE_DIR})

include(${CMAKE_SOURCE_DIR}/cmake/DownloadProject.cmake)

if (DOWNLOAD_GTEST)
    download_project(PROJ                googletest
            GIT_REPOSITORY      https://github.com/google/googletest.git
            GIT_TAG             master
            ${UPDATE_DISCONNECTED_IF_AVAILABLE}
            )

    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})

    # When using CMake 2.8.11 or later, header path dependencies
    # are automatically added to the gtest and gmock targets.
    # For earlier CMake versions, we have to explicitly add the
    # required directories to the header search path ourselves.
    if (CMAKE_VERSION VERSION_LESS 2.8.11)
        include_directories("${gtest_SOURCE_DIR}/include"
                "${gmock_SOURCE_DIR}/include")
    endif()
    set(GTEST_BOTH_LIBRARIES gtest gtest_main)
    set(GMOCK_BOTH_LIBRARIES gmock gmock_main)
else()
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
    message(${CMAKE_MODULE_PATH})
    find_package(GTest REQUIRED)
    find_package(GMock REQUIRED)
    include_directories(${GTEST_INCLUDE_DIRS})
    include_directories(${GMOCK_INCLUDE_DIRS})
endif()

add_custom_target(tests COMMENT "Builds all unit tests.")

function(mcl_add_test targetName targetFile)
    add_executable(${targetName} ${targetFile})
    target_link_libraries(${targetName} ${GTEST_BOTH_LIBRARIES} ${GMOCK_BOTH_LIBRARIES}
            mimiccomm ${CMAKE_THREAD_LIBS_INIT})
    add_test(NAME ${targetName}
            WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
            COMMAND ${targetName})
    add_dependencies(tests ${targetName})
endfunction()

function(mcl_add_mpi_test targetName targetFile numProcs)
    add_executable(${targetName} ${targetFile})
    target_link_libraries(${targetName} ${MPI_CXX_LIBRARIES}
            ${GTEST_BOTH_LIBRARIES} ${GMOCK_BOTH_LIBRARIES}
            mimiccomm ${CMAKE_THREAD_LIBS_INIT})
    set_target_properties(${targetName} PROPERTIES
            COMPILE_FLAGS "${MPI_CXX_COMPILER_FLAGS}"
            LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    add_test(NAME ${targetName}
            WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
            COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numProcs}
            ${MPIEXEC_PREFLAGS} ./${targetName} ${MPIEXEC_POSTFLAGS})
    add_dependencies(tests ${targetName})
endfunction()

mcl_add_test(api APITest.cpp)
mcl_add_test(endpoint EndpointsTest.cpp)
mcl_add_mpi_test(mpitransport MPITransportTest.cpp "3")