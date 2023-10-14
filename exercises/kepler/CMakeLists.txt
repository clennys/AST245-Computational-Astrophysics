cmake_minimum_required(VERSION 3.20)

project(
    kepler
    VERSION 1.0
    LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS yes)
set(CMAKE_C_STANDARD 23)
set(CMAKE_CXX_STANDARD 20)

find_package(fmt)

add_compile_options(
    # -Werror
    -Wall
    -Wextra
    -Wconversion
    -Wpedantic
    # note covered by above
    -Wshadow
    -Wundef
)

add_executable(kepler
    kepler.cpp
)

target_include_directories(kepler PUBLIC
    ../../includes/
)

target_link_libraries(kepler fmt::fmt)
