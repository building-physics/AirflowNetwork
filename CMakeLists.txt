cmake_minimum_required(VERSION 3.0.0)

project(afn)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(BUILD_TESTING "Build testing targets" ON)

if(CMAKE_COMPILER_IS_GNUCXX)
  # C++ 20 is allegedly available beginning in GCC v8, so some workarounds are needed
  # https://gcc.gnu.org/projects/cxx-status.html#cxx20
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "7.5.0")
    message(FATAL_ERROR "g++ versions earlier than 7.5.0 are not supported")
  endif()
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "10.0.0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDONT_HAVE_CPP20_NUMBERS")
    if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "10.0.0")
      set(BUILD_TESTING OFF)
    endif()
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
  endif()
elseif(MSVC)
  # https://devblogs.microsoft.com/cppblog/msvc-cpp20-and-the-std-cpp20-switch/
  if(${CMAKE_C_COMPILER_VERSION} VERSION_LESS "19.29.30040.0")
    message(FATAL_ERROR "Visual Studio earlier than VS2019 16.11 is not supported")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++20")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # clang build is untested
  # Need to check for clang version
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -stdlib=libc++")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/dependencies/include)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
