FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl
    GIT_TAG        cd73f3a4f3425fc70191f2c644c24de4642aa3a9
)
FetchContent_MakeAvailable(libigl)

add_executable(libigl_handle_deformation libigl_handle_deformation.cpp)
target_link_libraries(libigl_handle_deformation PRIVATE Mesh_smoothing_3 igl::core)

