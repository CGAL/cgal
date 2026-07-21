
FetchContent_Declare(
    cinolib
    GIT_REPOSITORY https://github.com/mlivesu/cinolib
    GIT_TAG        8fecf513ca7bf64556fcd73d7b7c2952778ad647
)
FetchContent_MakeAvailable(cinolib)

set(cinolib_DIR ${cinolib_SOURCE_DIR})
find_package(cinolib REQUIRED)


add_executable(cinolib_smooth cinolib_smooth.cpp)
target_link_libraries(cinolib_smooth PRIVATE Mesh_smoothing_3 cinolib)