cmake_minimum_required (VERSION 3.1...3.15)
project (DSM_REG)
set (targetNameA reg)
set (targetNameB transform)
set (targetNameC dsm_eval)

add_subdirectory(${targetNameA})
add_subdirectory(${targetNameB})
add_subdirectory(${targetNameC})

