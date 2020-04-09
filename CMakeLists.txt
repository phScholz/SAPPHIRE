#CMakeLists.txt

cmake_minimum_required (VERSION 3.1)
project (Sapphire)
message("Using CMAKE version ${CMAKE_VERSION})

if (GSL_PATH)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${GSL_PATH})
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${GSL_PATH})
endif (GSL_PATH) 

option (MPI_BUILD
       "Build Sapphire for MPI simulations. Note: Cross section calculations are disabled with MPI builds." OFF)

#Set the C++ standard to c++11
set (CMAKE_CXX_STANDARD 11)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

message("Checking if GNU SCIENTIFIC LIBRARY (GSL) can be found ...")
find_package (gsl REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

message("Checking if ROOT can be found ...")
find_package (ROOT REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})


if (MPI_BUILD)
add_definitions(-DMPI_BUILD)
message("Checking if Boost 1.4 COMPONENTS mpi serialization can be found ...")
FIND_PACKAGE(Boost 1.4 COMPONENTS mpi serialization REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
else (MPI_BUILD)
message("Checking if OpenMP can be found ...")
find_package (OpenMP REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (MPI_BUILD)

include_directories(${Sapphire_SOURCE_DIR}/coul/include)
add_subdirectory (coul)
set (EXTRA_LIBS coul)

configure_file( include/Sapphire_config.h.in ${CMAKE_BINARY_DIR}/generated/Sapphire_config.h )
include_directories( ${CMAKE_BINARY_DIR}/generated/ )

add_subdirectory (src)

