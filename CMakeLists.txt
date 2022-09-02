cmake_minimum_required(VERSION 3.12)

set(PROJ sirasac)
set(${PROJ}_VERSION "0.0.1")

project(${PROJ}
  VERSION "${${PROJ}_VERSION}"
  DESCRIPTION "Simple Random Sampling Clustering (SiRaSaC)."
  LANGUAGES C
)

add_subdirectory(c)
