cmake_minimum_required(VERSION 2.8)
project(TextureSynthesis)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

add_executable(TextureSynthesis test.cpp)
target_link_libraries(TextureSynthesis Aquila Ooura_fft)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    if(APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
endif()
