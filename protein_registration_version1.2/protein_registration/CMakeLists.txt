cmake_minimum_required(VERSION 3.9)
project(protein_registration)

set(CMAKE_CXX_STANDARD 17)

find_package( CGAL QUIET COMPONENTS  )
if ( NOT CGAL_FOUND )
    message(STATUS "This project requires the CGAL library, and will not be compiled.")
    return()
endif()

# include helper file
include( ${CGAL_USE_FILE} )
# Boost and its components
find_package( Boost REQUIRED COMPONENTS filesystem)
if ( NOT Boost_FOUND )
    message(STATUS "This project requires the Boost library, and will not be compiled.")
    return()
endif()


find_package (Eigen3 3.3 REQUIRED NO_MODULE)
# find_package(benchmark)
find_package(dlib)


# Setup your executable
# include_directories (BEFORE "../../include")
# include_directories (BEFORE "include")

set(CMAKE_CXX_FLAGS "-march=native")
set(CMAKE_CXX_FLAGS "-pthread")

add_executable(protein_registration src/main.cpp src/file_IO.h src/file_IO.cpp src/dSSD_sort.h src/DelaunayWrapper.h src/DelaunayWrapper.cpp src/hungarian_algorithm.h
        src/hungarian_algorithm.cpp src/SSD_sort.h src/ThreadPool.h src/configuration.cpp src/configuration.h)
#add_executable(protein_registration src/hungarian_algorithm.cpp)

target_link_libraries (protein_registration Eigen3::Eigen)
# target_link_libraries(protein_registration -lbenchmark)
target_link_libraries(protein_registration dlib::dlib)
target_link_libraries(protein_registration ${Boost_LIBRARIES})