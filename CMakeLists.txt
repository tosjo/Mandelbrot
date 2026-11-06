cmake_minimum_required(VERSION 3.16)

project(mandelbrot
    VERSION 2.0.0
    DESCRIPTION "Multithreaded Mandelbrot set renderer"
    LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(MANDELBROT_BUILD_TESTS "Build the unit tests" ON)
option(MANDELBROT_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

find_package(Threads REQUIRED)

add_library(mandelbrot_core STATIC
    src/options.cpp
    src/palette.cpp
    src/png.cpp
    src/render.cpp)
target_include_directories(mandelbrot_core PUBLIC include)
target_compile_features(mandelbrot_core PUBLIC cxx_std_17)
target_link_libraries(mandelbrot_core PUBLIC Threads::Threads)

add_library(mandelbrot_warnings INTERFACE)
if(MSVC)
    target_compile_options(mandelbrot_warnings INTERFACE /W4 /permissive-)
    if(MANDELBROT_WARNINGS_AS_ERRORS)
        target_compile_options(mandelbrot_warnings INTERFACE /WX)
    endif()
else()
    target_compile_options(mandelbrot_warnings INTERFACE
        -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow)
    if(MANDELBROT_WARNINGS_AS_ERRORS)
        target_compile_options(mandelbrot_warnings INTERFACE -Werror)
    endif()
endif()
target_link_libraries(mandelbrot_core PRIVATE mandelbrot_warnings)

add_executable(mandelbrot src/main.cpp)
target_compile_definitions(mandelbrot PRIVATE MANDELBROT_VERSION="${PROJECT_VERSION}")
target_link_libraries(mandelbrot PRIVATE mandelbrot_core mandelbrot_warnings)

if(MANDELBROT_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

include(GNUInstallDirs)
install(TARGETS mandelbrot RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
