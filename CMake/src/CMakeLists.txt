add_subdirectory_if( CGAL_USE_CGAL_CORE CGALCore)

add_subdirectory(CGAL)

optional_add_subdirectory( CGALimageIO OFF )
optional_add_subdirectory( CGALPDB     OFF )
optional_add_subdirectory( CGALQt      ON )

