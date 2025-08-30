cmake_minimum_required(VERSION 3.10.0)

project(ieee754_seq VERSION 0.1.0 LANGUAGES CXX)

project(ieee754_seq 
        VERSION 0.1.0
        DESCRIPTION "Header-only C++ class for defining ieee754 real number sequences"
        HOMEPAGE_URL https://gitlab.com/ir2-lab/ieee754_seq
        LANGUAGES CXX)
        
# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(GNUInstallDirs)

add_library(ieee754_seq INTERFACE)

target_include_directories(
  ieee754_seq 
  INTERFACE 
  ${CMAKE_INSTALL_INCLUDEDIR})

target_sources(ieee754_seq
  PUBLIC FILE_SET HEADERS
  BASE_DIRS ${PROJECT_SOURCE_DIR}/include
  FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/ieee754_seq.h
)

# Build test only if we are the main project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)

  add_executable(test EXCLUDE_FROM_ALL test/test.cpp )
  target_link_libraries(test ieee754_seq)

endif()

# install the headers that we specified with FILE_SET (header files for the user to include in their projects)
install(TARGETS ieee754_seq
  FILE_SET HEADERS
)