cmake_minimum_required(VERSION 3.20)

# set(VCPKG_TARGET_ARCHITECTURE arm64)
# set(VCPKG_CRT_LINKAGE dynamic)
# set(VCPKG_LIBRARY_LINKAGE dynamic)

set(CMAKE_TOOLCHAIN_FILE "/home/arminveres/gitfiles/vcpkg/scripts/buildsystems/vcpkg.cmake")

# set(CMAKE_C_COMPILER gcc-13)
# set(CMAKE_CXX_COMPILER g++-13)

project(n-body
    VERSION 0.1.0
    LANGUAGES C CXX
)

find_package(Eigen3 REQUIRED)
find_package(MathGL2 REQUIRED)
find_package(TBB REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS yes)
set(CMAKE_C_STANDARD 23)
set(CMAKE_CXX_STANDARD 20)

add_compile_options(
    -Wall
    -Wextra # reasonable and standard
    -Wshadow # warn the user if a variable declaration shadows one from a parent context
    -Wnon-virtual-dtor # warn the user if a class with virtual functions has a non-virtual destructor. This helps
    # catch hard to track down memory errors
    -Wold-style-cast # warn for c-style casts
    -Wcast-align # warn for potential performance problem casts
    -Wunused # warn on anything being unused
    -Woverloaded-virtual # warn if you overload (not override) a virtual function
    -Wpedantic # warn if non-standard C++ is used
    -Wconversion # warn on type conversions that may lose data
    -Wsign-conversion # warn on sign conversions
    -Wnull-dereference # warn if a null dereference is detected
    -Wdouble-promotion # warn if float is implicit promoted to double
    -Wformat=2 # warn on security issues around functions that format output (ie printf)
    -Wimplicit-fallthrough # warn on statements that fallthrough without an explicit annotation

    -Wmisleading-indentation # warn if indentation implies blocks where blocks do not exist
    -Wduplicated-cond # warn if if / else chain has duplicated conditions
    -Wduplicated-branches # warn if if / else branches have duplicated code
    -Wlogical-op # warn about logical operations being used where bitwise were probably wanted
    -Wuseless-cast # warn if you perform a cast to the same type

    # performance related flags
    # -mavx2
    -O3
    -march=native
    # -mfma
    # -mavx512f
    # -fopenmp

    # -O0
    # -g3
)

add_executable(${PROJECT_NAME}
    src/main.cpp
    src/particle.cpp
    src/data.cpp
    src/histogram.cpp
    src/shell.cpp
    src/system.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
    # "build/vcpkg_installed/arm64-osx/include/"
    "../../includes"
    "include"
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    $<IF:$<TARGET_EXISTS:mgl>,mgl,mgl-static>
    Eigen3::Eigen
    TBB::tbb
    # omp
)

