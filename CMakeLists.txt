cmake_minimum_required(VERSION 2.8)

PROJECT(sdoconv)

include(ExternalProject)

set(LOCAL_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/install )
set(LOCAL_PREFIX_PATH ${CMAKE_INSTALL_PREFIX} ${LOCAL_INSTALL_PREFIX} ${CMAKE_PREFIX_PATH})
string(REPLACE ";" "^^" LOCAL_PREFIX_PATH "${LOCAL_PREFIX_PATH}")

ExternalProject_Add (
  sdoconv
  LIST_SEPARATOR ^^
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src
  UPDATE_COMMAND ""
  INSTALL_COMMAND make install > ${LOCAL_INSTALL_PREFIX}/install_output.log
  INSTALL_DIR ${LOCAL_INSTALL_PREFIX}
  CMAKE_ARGS -DCMAKE_PREFIX_PATH=${LOCAL_PREFIX_PATH} -DCMAKE_INSTALL_PREFIX=${LOCAL_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  )

find_package( libsdo QUIET )
if( NOT libsdo_FOUND )
ExternalProject_Add(
  libsdo
  LIST_SEPARATOR ^^
  GIT_REPOSITORY https://github.com/rgottwald/libsdo.git
  GIT_TAG origin/master
  INSTALL_COMMAND ""
  CMAKE_ARGS -DCMAKE_PREFIX_PATH=${LOCAL_PREFIX_PATH} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
)
add_dependencies(sdoconv libsdo)
endif()

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

INSTALL( PROGRAMS ${LOCAL_INSTALL_PREFIX}/bin/sdoconv DESTINATION bin )
