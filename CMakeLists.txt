#CMakeLists.txt
message("\nStart building SAPPHIRE ... \n")

cmake_minimum_required(VERSION 3.1)
project(Sapphire C CXX)

set(CMAKE_CXX_STANDARD 11)

#Set optimization to 3
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

message("Using CMAKE version ${CMAKE_VERSION}")

set(build_types None Debug Regulus Release Travis MinSizeRel)

option(CODE_COVERAGE "Enable coverage reporting" OFF)
option(USE_GCOV "Create a GCov-enabled build." OFF)
option(BUILD_DOC "Build documentation" OFF)
option(GTEST "Build unit tests" ON)

message("CMAKE_CXX_COMPILER_ID     : ${CMAKE_CXX_COMPILER_ID}")
message("CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")


if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
  message("\nBuild type is set to ${CMAKE_BUILD_TYPE}\n")
  option(CODE_COVERAGE "Enable coverage reporting" OFF)
  message("OPTION CODE_COVERAGE IS OFF")
  option(USE_GCOV "Create a GCov-enabled build." OFF)
  message("OPTION USE_GCOV IS OFF")
  option(BUILD_DOC "Build documentation" OFF)
  message("OPTION BUILD_DOC IS OFF")
  option(GTEST "Build unit tests" ON)
  message("OPTION GTEST IS ON\n")
  enable_testing()  
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")

if(${CMAKE_BUILD_TYPE} MATCHES "Regulus")
  message("\nBuild type is set to ${CMAKE_BUILD_TYPE}\n")
  option(CODE_COVERAGE "Enable coverage reporting" OFF)
  message("OPTION CODE_COVERAGE IS OFF")
  option(USE_GCOV "Create a GCov-enabled build." OFF)
  message("OPTION USE_GCOV IS OFF")
  option(BUILD_DOC "Build documentation" OFF)
  message("OPTION BUILD_DOC IS OFF")
  option(GTEST "Build UNIT TESTS" OFF)
  message("OPTION GTEST IS OFF\n")  
endif(${CMAKE_BUILD_TYPE} MATCHES "Regulus")

if(${CMAKE_BUILD_TYPE} MATCHES "Travis")
  message("\nBuild type is set to ${CMAKE_BUILD_TYPE}\n")
  option(CODE_COVERAGE "Enable coverage reporting" ON)
  message("\nOPTION CODE_COVERAGE IS ON")
  option(USE_GCOV "Create a GCov-enabled build." ON)
  message("OPTION USE_GCOV IS ON")
  option(BUILD_DOC "Build documentation" ON)
  message("OPTION BUILD_DOC IS ON")
  option(GTEST "Build UNIT TESTS" ON)
  message("OPTION GTEST IS ON\n")
  enable_testing()  
endif(${CMAKE_BUILD_TYPE} MATCHES "Travis")


if (GSL_PATH)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${GSL_PATH})
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${GSL_PATH})
endif (GSL_PATH) 


#Set c++standard explicitly
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

#Explicitly setting link to math
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lm")

#Disable compilation warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w ")

# Code Coverage Configuration
add_library(coverage_config INTERFACE)

if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")

if (USE_GCOV)
  set(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")
  set(GCC_COVERAGE_LINK_FLAGS    "-lgcov")
endif()
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

message("\nChecking if GNU SCIENTIFIC LIBRARY (GSL) can be found ...")
find_package (gsl REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

message("\nChecking if ROOT can be found ...")
find_package (ROOT 6.0 REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})
message(STATUS "CXX_FLAGS                ${CMAKE_CXX_FLAGS}")

#add_definitions(-DMPI_BUILD)
message("\nChecking if Boost can be found ...")
FIND_PACKAGE(Boost 1.65 COMPONENTS filesystem random REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

message("\nChecking if OpenMP can be found ...")
find_package(OpenMP REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(${Sapphire_SOURCE_DIR}/coul/include)
add_subdirectory(coul)
set(EXTRA_LIBS coul)
configure_file(include/Sapphire_config.h.in ${CMAKE_BINARY_DIR}/generated/Sapphire_config.h )
include_directories( ${CMAKE_BINARY_DIR}/generated/ )

add_subdirectory(src)

if(${CMAKE_BUILD_TYPE} MATCHES "Travis" OR ${CMAKE_BUILD_TYPE} MATCHES "Release" OR ${CMAKE_BUILD_TYPE} MATCHES "Debug")
  add_subdirectory(tst)
endif(${CMAKE_BUILD_TYPE} MATCHES "Travis" OR ${CMAKE_BUILD_TYPE} MATCHES "Release" OR ${CMAKE_BUILD_TYPE} MATCHES "Debug")

  


if(BUILD_DOC MATCHES ON)
 # check if Doxygen is installed
 message("\nChecking if Doxygen can be found ...")
 find_package(Doxygen)
 if (DOXYGEN_FOUND)
     # set input and output files
     set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doxygen_config/sapphire.doxy)
     set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/sapphire.doxy)
 
     # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")
 
     # note the option ALL which allows to build the docs together with the application
      add_custom_target(doc_doxygen ALL
       COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
     COMMENT "Generating API documentation with Doxygen"
    VERBATIM )
 else (DOXYGEN_FOUND)
   message("Doxygen need to be installed to generate the doxygen documentation")
 endif (DOXYGEN_FOUND)
 endif(BUILD_DOC MATCHES ON)

message("Building Done\n")