cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Generate compile_commands.json for tools like clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compiler to clang
set(CMAKE_CXX_COMPILER clang++)

# Find all source files
file(GLOB_RECURSE SOURCES "src/*.cpp")
file(GLOB_RECURSE HEADERS "src/*.hpp")

# Add executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE src)

# Set compiler flags
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Werror)
