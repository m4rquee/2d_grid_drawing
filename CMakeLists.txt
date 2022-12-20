cmake_minimum_required(VERSION 3.23)
project(2d_draw)

set(CMAKE_CXX_STANDARD 23)

# Set a default build type if none was specified
set(default_build_type "Release")

# Headers directories
include_directories("include")
include_directories($ENV{GUROBI_HOME}/include)
# include_directories("lemon/lemon-1.3.1/include")

# Gurobi lib directory
link_directories($ENV{GUROBI_HOME}/lib)

# Lemon lib directory
# link_directories(${CMAKE_SOURCE_DIR}/lemon/lemon-1.3.1/lib)

# Compiler flags
# set(CMAKE_CXX_FLAGS "-Ofast -Wall -std=c++14 -D_GLIBCXX_USE_CXX11_ABI=0 -Wunused-local-typedefs -lemon")

# Target executable and linking
add_executable(2d_draw src/main.cpp)
target_link_libraries(2d_draw m gurobi_c++ gurobi95)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Running in debug mode")
    set_target_properties(2d_draw PROPERTIES COMPILE_FLAGS -DDEBUG_MODE)
endif()
