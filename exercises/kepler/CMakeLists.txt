cmake_minimum_required(VERSION 3.20)

project(
    kepler
    VERSION 1.0
    LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS yes)
set(CMAKE_C_STANDARD 23)
set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 REQUIRED)
find_package(Python3 REQUIRED)
find_package(fmt)

include_directories(SYSTEM
    ${PYTHON_INCLUDE_DIRS}
    "/usr/include/python3.12"
    "../../libs/matplotlib-cpp"
)

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

    -O3
    # -mavx2
    -mavx512f
    -march=native
)

add_executable(${PROJECT_NAME}
    main.cpp
    particle.cpp
    ode_solver.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
    "./includes"
    "../../includes"
    "../../libs/plotcpp"
)

target_link_libraries(${PROJECT_NAME}
    fmt::fmt
    Eigen3::Eigen
    ${PythonLibs_LIBRARIES}
    python3.12
)
