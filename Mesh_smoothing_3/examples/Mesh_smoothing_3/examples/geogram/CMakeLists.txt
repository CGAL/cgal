set(GEOGRAM_LIB_ONLY ON)

FetchContent_Declare(
    geogram
    GIT_REPOSITORY https://github.com/BrunoLevy/geogram
    GIT_TAG        72bfd854d4fb6bac360689a87fc7b10db87bb5f8 #master
)
FetchContent_MakeAvailable(geogram)


add_executable(geogram_smooth geogram_smooth.cpp)
add_executable(geogram_smooth_custom geogram_smooth_custom.cpp)
target_link_libraries(geogram_smooth PRIVATE Mesh_smoothing_3 geogram)
target_link_libraries(geogram_smooth_custom PRIVATE Mesh_smoothing_3 geogram)


