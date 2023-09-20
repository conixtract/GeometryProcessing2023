cmake_minimum_required(VERSION 3.8)
project(Assignment04)

add_executable(${PROJECT_NAME} "main.cc" "task.hh" "task.cc")

target_link_libraries(${PROJECT_NAME} PUBLIC
    typed-geometry
    polymesh
    glfw
    glow
    glow-extras
    ${COMMON_LINKER_FLAGS}
)
target_compile_options(${PROJECT_NAME} PUBLIC ${COMMON_COMPILER_FLAGS})
set_property(TARGET ${PROJECT_NAME} PROPERTY FOLDER "Assignment04")
