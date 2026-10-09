# Dépendances externes, compilées depuis les sources (sous-modules git, versions figées).
#   glfw   3.4              fenêtre + contexte OpenGL
#   glad   OpenGL 4.5 core  chargeur de fonctions (généré avec glad2, versionné ici)
#   imgui  1.92.9b-docking  interface (branche docking)
#   implot 1.0              graphes temps réel

# --- GLFW ---
set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL        OFF CACHE BOOL "" FORCE)
add_subdirectory(glfw)

# --- GLAD ---
add_library(glad STATIC glad/src/gl.c)
target_include_directories(glad SYSTEM PUBLIC glad/include)

# --- Dear ImGui ---
add_library(imgui STATIC
    imgui/imgui.cpp
    imgui/imgui_draw.cpp
    imgui/imgui_tables.cpp
    imgui/imgui_widgets.cpp
    imgui/imgui_demo.cpp
    imgui/backends/imgui_impl_glfw.cpp
    imgui/backends/imgui_impl_opengl3.cpp)
target_include_directories(imgui SYSTEM PUBLIC imgui imgui/backends)
target_link_libraries(imgui PUBLIC glfw)

# --- ImPlot ---
add_library(implot STATIC
    implot/implot.cpp
    implot/implot_items.cpp
    implot/implot_demo.cpp)
target_include_directories(implot SYSTEM PUBLIC implot)
target_link_libraries(implot PUBLIC imgui)
