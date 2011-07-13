include_directories(${PROJECT_SOURCE_DIR}/libs/f77)
btl_add_bench(btl_C main.cpp OFF)
# set_target_properties(btl_C PROPERTIES COMPILE_FLAGS "-fpeel-loops")