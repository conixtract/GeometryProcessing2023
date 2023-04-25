cmake_minimum_required(VERSION 3.8)
project(GeometryProcessingFramework)

# ===============================================
# global settings

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# ===============================================
# options

option(GP2_ENABLE_ASAN "if true, enables clang/MSVC address sanitizer" OFF)
option(GP2_ENABLE_MSAN "if true, enables clang/MSVC memory sanitizer" OFF)
option(GP2_ENABLE_UBSAN "if true, enables clang/MSVC undefined behaviour sanitizer" OFF)
option(GP2_ENABLE_TSAN "if true, enables clang/MSVC thread sanitizer" OFF)

if (GP2_ENABLE_ASAN AND GP2_ENABLE_TSAN)
    message(FATAL_ERROR "Can only enable one of TSan or ASan at a time")
endif()
if (GP2_ENABLE_ASAN AND GP2_ENABLE_MSAN)
    message(FATAL_ERROR "Can only enable one of ASan or MSan at a time")
endif()

option(GP2_ENABLE_WERROR "if true, enables -Werror, /WX" OFF)


# ===============================================
# compiler and linker flags

set(COMMON_COMPILER_FLAGS "")
set(COMMON_LINKER_FLAGS "")

if (MSVC)
    list(APPEND COMMON_COMPILER_FLAGS /MP)

    if (GP2_ENABLE_WERROR)
        list(APPEND COMMON_COMPILER_FLAGS /WX)
    endif()
else()
    list(APPEND COMMON_COMPILER_FLAGS -Wall -Wextra)

    if (GP2_ENABLE_WERROR)
        list(APPEND COMMON_COMPILER_FLAGS -Werror)
    endif()

    if (GP2_ENABLE_ASAN OR GP2_ENABLE_TSAN OR GP2_ENABLE_MSAN OR GP2_ENABLE_UBSAN)
        list(APPEND COMMON_COMPILER_FLAGS -fno-omit-frame-pointer -g)
        list(APPEND COMMON_LINKER_FLAGS -fno-omit-frame-pointer -g)
    endif()

    if (GP2_ENABLE_ASAN)
        list(APPEND COMMON_COMPILER_FLAGS -fsanitize=address)
        list(APPEND COMMON_LINKER_FLAGS -fsanitize=address)
    endif()

    if (GP2_ENABLE_TSAN)
        list(APPEND COMMON_COMPILER_FLAGS -fsanitize=thread)
        list(APPEND COMMON_LINKER_FLAGS -fsanitize=thread)
    endif()

    if (GP2_ENABLE_MSAN)
        list(APPEND COMMON_COMPILER_FLAGS -fsanitize=memory)
        list(APPEND COMMON_LINKER_FLAGS -fsanitize=memory)
    endif()

    if (GP2_ENABLE_UBSAN)
        list(APPEND COMMON_COMPILER_FLAGS
            -fsanitize=undefined
            -fno-sanitize-recover=all
            -fno-sanitize=alignment,vptr
        )
        list(APPEND COMMON_LINKER_FLAGS
            -fsanitize=undefined
            -fno-sanitize-recover=all
            -fno-sanitize=alignment,vptr
        )
    endif()
endif()

# ===============================================
# Bin dir
if(MSVC)
    set(BIN_DIR ${CMAKE_SOURCE_DIR}/bin)
elseif(CMAKE_BUILD_TYPE STREQUAL "")
    set(BIN_DIR ${CMAKE_SOURCE_DIR}/bin/Default)
else()
    set(BIN_DIR ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${BIN_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${BIN_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${BIN_DIR})
set(GLOW_BIN_DIR ${CMAKE_SOURCE_DIR}/bin)

# ===============================================
# disable glfw additionals
option(GLFW_BUILD_EXAMPLES  OFF)
option(GLFW_BUILD_TESTS  OFF)
option(GLFW_BUILD_DOCS  OFF)
option(GLFW_INSTALL  OFF)


# ===============================================s
# add submodules
add_subdirectory(extern/clean-core)
add_subdirectory(extern/typed-geometry)
add_subdirectory(extern/polymesh)
add_subdirectory(extern/glfw)
add_subdirectory(extern/glow)
add_subdirectory(extern/imgui)
add_subdirectory(extern/glow-extras)
add_subdirectory(extern/eigen-lean)

# ===============================================
# configure Source grouping
file(GLOB_RECURSE SOURCES
    "src/*.cc"
    "src/*.hh"
    "src/*.cpp"
    "src/*.h"
)
# group sources according to folder structure
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCES})


# ===============================================
# dependency grouping in the IDE
foreach(TARGET_NAME
    clean-core
    typed-geometry
    polymesh
    glfw
    glow
    glow-extras
    )
    set_property(TARGET ${TARGET_NAME} PROPERTY FOLDER "Extern")
endforeach()


# ===============================================
# add executable

file(GLOB_RECURSE mains "src/*/main.cc")
foreach(main ${mains})
    get_filename_component(path ${main} DIRECTORY)
    string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" spath ${path})
    message(STATUS "adding exercise ${spath}")
    add_subdirectory(${path})
endforeach()


