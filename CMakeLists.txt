cmake_minimum_required(VERSION 3.15)
project(mces VERSION 1.0.0 LANGUAGES C)

# Configure C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Common compiler flags - defined once, used everywhere
set(COMMON_COMPILE_FLAGS
    $<$<CONFIG:Debug>:-Wall -Wextra -g>
    $<$<CONFIG:Release>:-Wall -Wextra -O3 -DNDEBUG>
)

add_executable(mces src/main.c src/code_source_1d.c)
target_include_directories(mces PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(mces PRIVATE m)
