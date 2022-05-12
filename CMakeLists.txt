cmake_minimum_required(VERSION 3.21)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/.cmake/modules/")

include(CheckIPOSupported)

project(samtag)

set(samtag_VERSION_MAJOR 1)
set(samtag_VERSION_MINOR 0)

# Generate list of compile commands.  This helps debugging and doesn't have a downside.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

configure_file (
    "${PROJECT_SOURCE_DIR}/src/version.h.in"
    "${PROJECT_BINARY_DIR}/generated/version.hpp"
)

include_directories(${CMAKE_BINARY_DIR}/generated)

set(COMPILER_ARCHITECTURE "native" CACHE STRING "Compiler -march argument")

set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message("-- Build type: " ${CMAKE_BUILD_TYPE})

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${samtag_SOURCE_DIR}" CACHE PATH "Installing to ${samtag_SOURCE_DIR}$" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
message("-- Installation prefix: " ${CMAKE_INSTALL_PREFIX})

add_subdirectory(src)
