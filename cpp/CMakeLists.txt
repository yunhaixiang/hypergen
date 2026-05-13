cmake_minimum_required(VERSION 3.20)
project(hyperelliptic_data_gen_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(PkgConfig REQUIRED)
pkg_check_modules(FLINT REQUIRED flint)
pkg_check_modules(SQLITE3 REQUIRED sqlite3)

add_executable(hyperelliptic_cpp hyperelliptic.cpp)
target_include_directories(hyperelliptic_cpp PRIVATE ${FLINT_INCLUDE_DIRS} ${SQLITE3_INCLUDE_DIRS})
target_link_directories(hyperelliptic_cpp PRIVATE ${FLINT_LIBRARY_DIRS} ${SQLITE3_LIBRARY_DIRS})
target_link_libraries(hyperelliptic_cpp PRIVATE ${FLINT_LIBRARIES} ${SQLITE3_LIBRARIES})
target_compile_options(hyperelliptic_cpp PRIVATE ${FLINT_CFLAGS_OTHER} ${SQLITE3_CFLAGS_OTHER} -Wall -Wextra -Wpedantic)
