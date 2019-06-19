# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Polygonal_surface_reconstruction_Examples )


cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

include_directories( BEFORE ../../include )

# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "polyfit_example_without_input_planes.cpp" )
create_single_source_cgal_program( "polyfit_example_user_provided_planes.cpp" )
create_single_source_cgal_program( "polyfit_example_model_complexty_control.cpp" )


find_package( SCIP )

if ( NOT SCIP_FOUND )
    find_package( GLPK )

    if ( NOT GLPK_FOUND )
      message(STATUS "This project requires either SCIP or GLPK, and will not be compiled. "
                      "Please proivde either 'SCIP_DIR' or 'GLPK_INCLUDE_DIR' and 'GLPK_LIBRARIES'")
    else()
        add_compile_definitions(CGAL_USE_GLPK)

        include_directories( BEFORE ${GLPK_INCLUDE_DIR} )

        target_link_libraries( polyfit_example_without_input_planes PRIVATE ${GLPK_LIBRARIES} )
        target_link_libraries( polyfit_example_user_provided_planes PRIVATE ${GLPK_LIBRARIES} )
        target_link_libraries( polyfit_example_model_complexty_control PRIVATE ${GLPK_LIBRARIES} )

        message("GLPK found and used")

    endif()

else()

    add_compile_definitions(CGAL_USE_SCIP)

    include_directories( BEFORE ${SCIP_INCLUDE_DIRS} )

    target_link_libraries( polyfit_example_without_input_planes PRIVATE ${SCIP_LIBRARIES} )
    target_link_libraries( polyfit_example_user_provided_planes PRIVATE ${SCIP_LIBRARIES} )
    target_link_libraries( polyfit_example_model_complexty_control PRIVATE ${SCIP_LIBRARIES} )

    message("SCIP found and used")

endif()

