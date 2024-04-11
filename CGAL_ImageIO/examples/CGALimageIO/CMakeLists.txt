# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(CGALimageIO_Examples)

find_package(CGAL REQUIRED COMPONENTS ImageIO)

create_single_source_cgal_program("convert_raw_image_to_inr.cpp")
create_single_source_cgal_program("test_imageio.cpp")
create_single_source_cgal_program("extract_a_sub_image.cpp")
create_single_source_cgal_program("slice_image.cpp")
