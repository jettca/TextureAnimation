cmake_minimum_required(VERSION 2.8)
project(TextureSynthesis)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

add_executable(TextureSynthesis
    src/test.cpp
    src/Synthesis/AudioDevice.cpp src/Synthesis/AudioDevice.h
    src/Synthesis/StatsGenerator.h src/Synthesis/StatsGenerator.cpp
    src/Synthesis/TextureFilterer.h src/Synthesis/TextureFilterer.cpp
    src/Synthesis/TextureSynthesizer.h src/Synthesis/TextureSynthesizer.cpp
    src/Filtering/FilterBank.cpp src/Filtering/FilterBank.h
    src/Filtering/Filter.cpp src/Filtering/Filter.h
    src/Filtering/LowpassFilter.cpp src/Filtering/LowpassFilter.h
    src/Filtering/CochlearFilter.cpp src/Filtering/CochlearFilter.h
    src/Filtering/Signal.h src/Filtering/Signal.cpp
    src/Filtering/Downsampler.h src/Filtering/Downsampler.cpp
)

include(FindPkgConfig)
pkg_search_module(SDL2 REQUIRED sdl2)
pkg_search_module(REQUIRED SDL2_mixer REQUIRED)

target_link_libraries(TextureSynthesis
    Aquila
    Ooura_fft
    gsl
    gslcblas
    ${SDL2_LIBRARIES}
)

include_directories(
    src
)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    if(APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif()
endif()
