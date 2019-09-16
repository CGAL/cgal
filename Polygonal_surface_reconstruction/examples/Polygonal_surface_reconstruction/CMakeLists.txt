# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Polygonal_surface_reconstruction_Examples )


cmake_minimum_required(VERSION 3.1...3.15)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "NOTICE: This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "NOTICE: This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "polyfit_example_without_input_planes.cpp" )
create_single_source_cgal_program( "polyfit_example_user_provided_planes.cpp" )
create_single_source_cgal_program( "polyfit_example_model_complexty_control.cpp" )
create_single_source_cgal_program( "polyfit_example_with_region_growing.cpp" )

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
if(EIGEN3_FOUND)
  include( ${EIGEN3_USE_FILE} )

  find_package( SCIP QUIET)

  if ( NOT SCIP_FOUND )
      find_package( GLPK QUIET)

      if ( NOT GLPK_FOUND )
        message( STATUS "NOTICE: This project requires either SCIP or GLPK, and will not be compiled. "
                        "Please provide either 'SCIP_DIR' or 'GLPK_INCLUDE_DIR' and 'GLPK_LIBRARIES'")
      else()
        

          include_directories( BEFORE ${GLPK_INCLUDE_DIR} )

          target_link_libraries( polyfit_example_without_input_planes PRIVATE ${GLPK_LIBRARIES} )
          target_compile_definitions(polyfit_example_without_input_planes PRIVATE -DCGAL_USE_GLPK)
          target_link_libraries( polyfit_example_user_provided_planes PRIVATE ${GLPK_LIBRARIES} )
          target_compile_definitions(polyfit_example_user_provided_planes PRIVATE -DCGAL_USE_GLPK)
          target_link_libraries( polyfit_example_model_complexty_control PRIVATE ${GLPK_LIBRARIES} )
          target_compile_definitions(polyfit_example_model_complexty_control PRIVATE -DCGAL_USE_GLPK)
          target_link_libraries( polyfit_example_with_region_growing PRIVATE ${GLPK_LIBRARIES} )
          target_compile_definitions(polyfit_example_with_region_growing PRIVATE -DCGAL_USE_GLPK)

          message("GLPK found and used")

      endif()

  else()

      include_directories( BEFORE ${SCIP_INCLUDE_DIRS} )

      target_link_libraries( polyfit_example_without_input_planes PRIVATE ${SCIP_LIBRARIES} )
      target_compile_definitions(polyfit_example_without_input_planes PRIVATE -DCGAL_USE_SCIP)
      target_link_libraries( polyfit_example_user_provided_planes PRIVATE ${SCIP_LIBRARIES} )
      target_compile_definitions(polyfit_example_user_provided_planes PRIVATE -DCGAL_USE_SCIP)
      target_link_libraries( polyfit_example_model_complexty_control PRIVATE ${SCIP_LIBRARIES} )
      target_compile_definitions(polyfit_example_model_complexty_control PRIVATE -DCGAL_USE_SCIP)
      target_link_libraries( polyfit_example_with_region_growing PRIVATE ${SCIP_LIBRARIES} )
      target_compile_definitions(polyfit_example_with_region_growing PRIVATE -DCGAL_USE_SCIP)

      message("SCIP found and used")

  endif()

else()
  message(STATUS "NOTICE: Some of the executables in this directory need Eigen 3.1 (or greater) and will not be compiled.")
endif()
