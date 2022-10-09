cmake_minimum_required(VERSION 3.2)

project(torchtest)

set(CMAKE_PREFIX_PATH /mnt/c/Users/ENCHOU/Documents/git-repo/Fractional-Step-FDM-Staggered-Lid-Driven-Cavity-/libtorch)

find_package(Torch REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} "${TORCH_LIBRARIES}")
target_link_libraries(${PROJECT_NAME} "${TORCH_CUDA_LIBRARIES}")

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 14)

