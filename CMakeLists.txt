cmake_minimum_required(VERSION 3.20)
project(PhysicsLab3D VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Type de build" FORCE)
endif()

option(PHYSICSLAB_BUILD_TESTS "Compiler les tests de validation" ON)

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/third_party/glfw/CMakeLists.txt")
    message(FATAL_ERROR "Sous-modules absents. Lancer : git submodule update --init --recursive")
endif()

add_subdirectory(third_party)

# Avertissements stricts pour notre code uniquement (pas pour third_party).
function(physicslab_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endfunction()

# --- Cœur physique : maths, constantes, solveurs, monde (aucune dépendance graphique) ---
add_library(physicslab_core STATIC
    src/core/Solver.cpp
    src/core/Events.cpp
    src/core/World.cpp
    src/mechanics/Projectile.cpp
    src/mechanics/Oscillator.cpp
    src/mechanics/Pendulum.cpp
    src/mechanics/DoublePendulum.cpp
    src/mechanics/Kepler.cpp
    src/mechanics/NBody.cpp
    src/mechanics/Friction.cpp
    src/mechanics/Collision.cpp
    src/mechanics/Bounce.cpp
    src/mechanics/Cradle.cpp)
target_include_directories(physicslab_core PUBLIC include)
physicslab_warnings(physicslab_core)

# --- Rendu OpenGL (découplé de la physique) ---
add_library(physicslab_render STATIC
    src/render/Camera.cpp
    src/render/Renderer.cpp)
target_include_directories(physicslab_render PUBLIC include)
target_link_libraries(physicslab_render PUBLIC glad)
physicslab_warnings(physicslab_render)

# --- Application ---
add_executable(physicslab
    src/app/main.cpp
    src/app/Application.cpp
    src/app/UiCommon.cpp
    src/app/ProjectileModule.cpp
    src/app/OscillatorModule.cpp
    src/app/PendulumModule.cpp
    src/app/DoublePendulumModule.cpp
    src/app/KeplerModule.cpp
    src/app/NBodyModule.cpp
    src/app/FrictionModule.cpp
    src/app/BounceModule.cpp
    src/app/CradleModule.cpp)
target_link_libraries(physicslab PRIVATE physicslab_core physicslab_render glfw imgui implot)
target_compile_definitions(physicslab PRIVATE
    PHYSICSLAB_SHADER_DIR="${CMAKE_SOURCE_DIR}/shaders")
physicslab_warnings(physicslab)

# --- Tests ---
if(PHYSICSLAB_BUILD_TESTS)
    enable_testing()
    add_executable(test_core tests/test_core.cpp)
    target_link_libraries(test_core PRIVATE physicslab_core)
    physicslab_warnings(test_core)
    add_test(NAME core COMMAND test_core)
endif()
