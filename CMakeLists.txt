cmake_minimum_required(VERSION 3.10)
project(TASK)

set(CMAKE_CXX_STANDARD 14)

set(UTIL_SRC util/trie.cc util/util.cc util/PB.cpp util/ioUint.cpp util/ioPB.cpp util/mtrie.cc util/canLabel.cpp util/canLabelET.cpp)

set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -g -Wall -Wno-unused-variable -pthread")
set(CMAKE_CXX_FLAGS_RELEASE   "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -O3 -pthread -fopenmp")

add_executable(buildR2T buildR2T.cpp ${UTIL_SRC})
add_executable(IUA IUA.cpp ${UTIL_SRC})
add_executable(IQA IQA.cpp ${UTIL_SRC})
add_executable(multiKey multiKey.cpp ${UTIL_SRC})

FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
