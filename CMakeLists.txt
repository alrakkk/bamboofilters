cmake_minimum_required(VERSION 3.5)

project(bamboofilters)

set(CMAKE_CXX_FLAGS
    "-fno-strict-aliasing -O3 -std=c++11 -march=core-avx2 ${CMAKE_CXX_FLAGS} -fopenmp"
)
set(CMAKE_EXE_LINKER_FLAGS "-lpthread ${CMAKE_EXE_LINKER_FLAGS}")
set(CMAKE_CXX_FLAGS "-fpermissive ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_subdirectory(src)
add_subdirectory(test)

add_subdirectory(tools)

# Add download target
add_custom_target(download_ecoli
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/download_ecoli.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Downloading E. coli reference genome"
)
