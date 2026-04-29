set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
cmake_minimum_required (VERSION 3.1...3.15)
set (targetName dsm_eval)

add_executable (${targetName} dsm_eval.cpp ../src/ICP.cpp ../src/kdtree.cpp ../src/example.cpp ../src/svd.cpp ../src/argparse.hpp)
include_directories( ../src)
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(${targetName} PRIVATE Eigen3::Eigen)

#GDAL
find_package(GDAL REQUIRED)
target_include_directories(${targetName} PRIVATE ${GDAL_INCLUDE_DIRS})
target_link_libraries(${targetName} PRIVATE ${GDAL_LIBRARIES})
