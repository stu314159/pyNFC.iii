cmake_minimum_required (VERSION 3.10)
project(pyNFCiii VERSION 0.1 LANGUAGES CXX)

find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(PythonLibs 2.7 REQUIRED)
find_package(BOOST COMPONENTS python REQUIRED)

add_library(PartitionHelper MODULE PartitionHelper.cpp)

target_link_libraries(PartitionHelper ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
target_include_directories(PartitionHelper PUBLIC include
                                           PRIVATE ${PYTHON_INCLUDE_DIRS})
