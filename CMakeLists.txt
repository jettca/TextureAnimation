cmake_minimum_required(VERSION 2.8)
project(TextureSynthesis)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

add_executable(TextureSynthesis
    src/test.cpp
    src/AudioDevice.cpp src/AudioDevice.h
    src/FilterBank.cpp src/FilterBank.h
    src/Filter.cpp src/Filter.h
    src/LowpassFilter.cpp src/LowpassFilter.h
    src/CochlearFilter.cpp src/CochlearFilter.h
    src/Signal.h src/Signal.cpp
    src/AuditoryModel.h src/AuditoryModel.cpp
)


include(FindPkgConfig)
pkg_search_module(SDL2 REQUIRED sdl2)
pkg_search_module(REQUIRED SDL2_mixer REQUIRED)

target_link_libraries(TextureSynthesis
    Aquila
    Ooura_fft
    ${SDL2_LIBRARIES}
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
