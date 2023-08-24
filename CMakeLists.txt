
cmake_minimum_required(VERSION 3.12)

project(makeSEAL VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_BUILD_TYPE debug)

# Enable OpenMP if available
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Enable if doing valgrind memory check
# include(CTest)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g ")
# enable_testing()

add_subdirectory(3rdparty)
add_subdirectory(run)