# MIT License
# 
# Copyright (c) 2022 Dr. Thomas Lang
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

set(SIRASAC_C sirasac)

# ----------------------------------------------------------------------------------------- #
option(${SIRASAC_C}_SHARED     "Build shared library?"          ON)
option(${SIRASAC_C}_BUILD_APP  "Build standalone application?" OFF)
option(SIRASAC_PYTHON_BINDINGS "Build Python bindings."         ON)

if(${SIRASAC_C}_SHARED)
  set(${SIRASAC_C}_BUILD_TYPE SHARED)
else()
  set(${SIRASAC_C}_BUILD_TYPE STATIC)
endif()

if(SIRASAC_PYTHON_BINDINGS)
  set(${SIRASAC_C}_BUILD_TYPE STATIC CACHE BOOL "Building Python library; Building ${SIRASAC_C} as static version." FORCE)
endif()

message(STATUS "[${SIRASAC_C}] Building library in ${${SIRASAC_C}_BUILD_TYPE} version.")

# ----------------------------------------------------------------------------------------- #

set(src ${CMAKE_CURRENT_SOURCE_DIR})
set(SAMPLING_SOURCES
  ${src}/RandomSampling.h
  ${src}/RandomSampling.c
  ${src}/RandomNumber.h
  ${src}/RandomNumber.c
)
set(EXPORT_SOURCES
  ${src}/sirasac.h
  ${src}/sirasac.c
)
set(_ALL_LIBRARY_SOURCES
  ${SAMPLING_SOURCES}
  ${EXPORT_SOURCES}
)
set(EXPORTED_HEADERS
  ${src}/sirasac.h
)

# ----------------------------------------------------------------------------------------- #
# Create the library.
add_library(${SIRASAC_C} ${${SIRASAC_C}_BUILD_TYPE} ${_ALL_LIBRARY_SOURCES})

set_target_properties(${SIRASAC_C} PROPERTIES PUBLIC_HEADER "${EXPORTED_HEADERS}")
set_target_properties(${SIRASAC_C} PROPERTIES VERSION ${${PROJ}_VERSION})
set_target_properties(${SIRASAC_C} PROPERTIES LINKER_LANGUAGE C)
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Build standalone application if requested.
if(${SIRASAC_C}_BUILD_APP)
  add_executable(${SIRASAC_C}_APP
    ${src}/app/main.c
    ${EXPORTED_HEADERS}
  )
  target_link_libraries(${SIRASAC_C}_APP PUBLIC ${SIRASAC_C})
endif()
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Provide installation rule, upon invocation this will install the header files into 
# ${CMAKE_INSTALL_PREFIX}/include  and the library into  ${CMAKE_INSTALL_PREFIX}/lib
install(TARGETS ${SIRASAC_C}
  LIBRARY       DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
  ARCHIVE       DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
  RUNTIME       DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
)
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Find OpenMP and use it.
find_package(OpenMP)
if(OpenMP_C_FOUND)
  target_link_libraries(${SIRASAC_C} PUBLIC OpenMP::OpenMP_C)

  if(${SIRASAC_C}_BUILD_APP)
    target_link_libraries(${SIRASAC_C}_APP PUBLIC OpenMP::OpenMP_C)
  endif()
endif()
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Generate the export header for the library.
include(GenerateExportHeader)
generate_export_header(${SIRASAC_C}
  BASE_NAME         sirasac
  EXPORT_MACRO_NAME sirasac_API
  EXPORT_FILE_NAME  ${CMAKE_INSTALL_PREFIX}/include/sirasac_export.h
)
target_include_directories(${SIRASAC_C} PUBLIC ${CMAKE_INSTALL_PREFIX}/include)
# ----------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------- #
# Create nice mappings for IDEs.
source_group("sirasac"           FILES ${EXPORT_SOURCES})
source_group("sirasac\\Sampling" FILES ${SAMPLING_SOURCES})

if(${SIRASAC_C}_BUILD_APP)
  source_group("sirasac\\Application" FILES ${src}/app/main.c)
endif()
# ----------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------- #
# Optional: Build Python bindings for C library as well.

if(SIRASAC_PYTHON_BINDINGS)
  find_package(PythonInterp 3.7 REQUIRED)
  find_package(PythonLibs 3.7 REQUIRED)
  
  # Variables to be substituted in the Python binding setup.py file.
  set(SIRASAC_C_HEADER_DIR        "${src}")
  set(SIRASAC_C_LIBRARY_DIR       "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>")
  set(SIRASAC_C_EXPORT_HEADER_DIR "${CMAKE_INSTALL_PREFIX}/include")
  set(SIRASAC_PY_SOURCE_DIR       "${CMAKE_CURRENT_SOURCE_DIR}/../python")
  set(SIRASAC_VERSION             ${${PROJ}_VERSION})
  set(SIRASAC_PY_INSTALL_DIR      ${CMAKE_INSTALL_PREFIX})
  
  set(SIRASAC_PY_SETUP "${CMAKE_CURRENT_BINARY_DIR}/setup_$<CONFIG>.py")
  set(SIRASAC_PY_SETUP_TEMPLATE "${SIRASAC_PY_SOURCE_DIR}/setup.py.in")
  
  file(READ ${SIRASAC_PY_SETUP_TEMPLATE} setupContents)
  string(CONFIGURE "${setupContents}" setupContents @ONLY)
  file(GENERATE
    OUTPUT ${SIRASAC_PY_SETUP}
    CONTENT "${setupContents}"
  )

  add_custom_command(
    TARGET ${SIRASAC_C} POST_BUILD
    COMMAND ${PYTHON_EXECUTABLE} ${SIRASAC_PY_SETUP} bdist_wheel -d ${SIRASAC_PY_INSTALL_DIR}
    COMMENT "Building Python bindings."
    DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/../python/c_random_sampling.pxd
      ${CMAKE_CURRENT_SOURCE_DIR}/../python/random_sampling.pyx
  )

endif()

# ----------------------------------------------------------------------------------------- #
