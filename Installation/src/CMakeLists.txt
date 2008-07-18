add_subdirectory_if( CGAL_USE_CGAL_CORE CGALCore)

add_subdirectory(CGAL)

optional_add_subdirectory( CGALimageIO ON )
optional_add_subdirectory( CGALPDB     ON )
optional_add_subdirectory( CGALQt      ON )

