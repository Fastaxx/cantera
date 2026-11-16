# Builds the shared Cantera library, including its C API (CLib), from the source
# distribution assembled by build_sdist.py. The CLib sources and the Julia bindings
# are already generated, so only CMake, a C++20 compiler and the third-party
# libraries are needed; this is what the Yggdrasil recipe for Cantera_jll builds.
cmake_minimum_required(VERSION 3.21...3.31)
project(cantera_clib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(yaml-cpp REQUIRED)
find_package(SUNDIALS REQUIRED)
find_package(fmt REQUIRED)
# Only Boost's header-only parts are used. Some Boost installations (e.g. Yggdrasil's
# for Windows) ship the headers without CMake package files, so fall back to them.
find_package(Boost 1.70 CONFIG QUIET)
if(NOT Boost_FOUND)
  find_path(BOOST_INCLUDE_DIR boost/version.hpp REQUIRED)
  add_library(Boost::headers INTERFACE IMPORTED)
  set_target_properties(Boost::headers PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${BOOST_INCLUDE_DIR}")
endif()

# Optional features of the SCons build that this library does not use.
set(CT_USE_LAPACK 0)
set(CT_SUNDIALS_USE_LAPACK 0)
set(CT_USE_HDF5 0)
set(CT_USE_SYSTEM_HIGHFIVE 0)
configure_file(include/cantera/base/config.h.in
               "${CMAKE_CURRENT_BINARY_DIR}/include/cantera/base/config.h")

file(GLOB_RECURSE CT_SOURCES CONFIGURE_DEPENDS src/*.cpp clib/src/*.cpp)
add_library(cantera SHARED ${CT_SOURCES})
# ct_defs.h includes "config.h" by its bare name, hence the extra directory.
target_include_directories(cantera PRIVATE
  "${CMAKE_CURRENT_BINARY_DIR}/include" "${CMAKE_CURRENT_BINARY_DIR}/include/cantera/base"
  include src src/base clib/include)
target_link_libraries(cantera PRIVATE
  yaml-cpp::yaml-cpp Eigen3::Eigen fmt::fmt Boost::headers
  SUNDIALS::core SUNDIALS::nvecserial SUNDIALS::cvodes SUNDIALS::idas
  SUNDIALS::sunmatrixband SUNDIALS::sunmatrixdense SUNDIALS::sunmatrixsparse
  SUNDIALS::sunlinsolband SUNDIALS::sunlinsoldense SUNDIALS::sunlinsolspgmr
  SUNDIALS::sunnonlinsolnewton)
set_target_properties(cantera PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  # The headers do not annotate the API with dllexport; export every symbol.
  WINDOWS_EXPORT_ALL_SYMBOLS ON)

include(GNUInstallDirs)
install(TARGETS cantera
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(DIRECTORY clib/include/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
install(DIRECTORY data/ DESTINATION "${CMAKE_INSTALL_DATADIR}/cantera/data")
install(DIRECTORY julia/ DESTINATION "${CMAKE_INSTALL_DATADIR}/cantera/julia")
