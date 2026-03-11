cmake_minimum_required(VERSION 3.24)

project(voronoi-partition VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(cmake/CPM.cmake)

CPMAddPackage("gh:CLIUtils/CLI11#v2.6.1")
CPMAddPackage("gh:liukaizheng/gpf-algorithm#main")
find_package(CGAL REQUIRED)
if(NOT TARGET CGAL::CGAL_Core)
    add_library(CGAL::CGAL_Core INTERFACE IMPORTED)
    set_target_properties(CGAL::CGAL_Core PROPERTIES
        INTERFACE_LINK_LIBRARIES CGAL::CGAL
    )
endif()
find_package(Eigen3 REQUIRED)
CPMAddPackage(
    NAME libigl
    GIT_TAG main
    GITHUB_REPOSITORY liukaizheng/libigl
    OPTIONS
        "LIBIGL_COPYLEFT_CGAL ON"
)

add_executable(voronoi-partition
    src/main.cpp
)

target_compile_features(voronoi-partition PRIVATE cxx_std_23)
target_link_libraries(voronoi-partition PRIVATE
    CLI11::CLI11
    gpf::algorithm
    CGAL::CGAL
    igl::core
    igl_copyleft::cgal
)
target_include_directories(voronoi-partition PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_options(voronoi-partition PRIVATE
    # -Wall -Wextra -Wpedantic
    # -fsanitize=address -fno-omit-frame-pointer
)

target_link_options(voronoi-partition PRIVATE
    # -fsanitize=address
)
