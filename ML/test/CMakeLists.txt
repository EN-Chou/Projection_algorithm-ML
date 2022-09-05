cmake_minimum_required(VERSION 3.2)

project(torchtest)

<<<<<<< HEAD
set(CMAKE_PREFIX_PATH /mnt/c/Users/ENCHOU/Documents/git-repo/Fractional-Step-FDM-Staggered-Lid-Driven-Cavity-/libtorch)
=======
set(CMAKE_PREFIX_PATH /home/enchou/Downloads/libtorch)
>>>>>>> 7d55f6f9ed629caa6f5fd8cf7a29fbac6a2d626e

find_package(Torch REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

add_executable(${PROJECT_NAME} test.cpp)

target_link_libraries(${PROJECT_NAME} "${TORCH_LIBRARIES}")

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 14)

