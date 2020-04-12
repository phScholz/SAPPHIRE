#CMakeLists.txt

cmake_minimum_required (VERSION 3.1)
project (Sapphire)
message("Using CMAKE version ${CMAKE_VERSION}")

if (GSL_PATH)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${GSL_PATH})
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${GSL_PATH})
endif (GSL_PATH) 

option (MPI_BUILD "Build Sapphire for MPI simulations. Note: Cross section calculations are disabled with MPI builds." OFF)

#Set the C++ standard to c++11
set (CMAKE_CXX_STANDARD 11)

#Set optimization to 3
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -O3")

#Set compilation warnings
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

message("Checking if GNU SCIENTIFIC LIBRARY (GSL) can be found ...")
find_package (gsl REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

message("Checking if ROOT can be found ...")
find_package (ROOT REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

set(CMAKE_BUILD_PARALLEL_LEVEL "4")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -j 20")

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

## first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

# check if Doxygen is installed
message("Checking if Doxygen can be found ...")
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doxygen_config/sapphire.doxy)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
   configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
   message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
  add_custom_target( doc_doxygen ALL
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
   VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

message("Done")
