message("Configuring libCGAL_Core")

collect_cgal_library(CGAL_Core "")

if(CGAL_HEADER_ONLY)
  set(keyword INTERFACE)
else()
  set(keyword PUBLIC)
endif()

target_link_libraries( CGAL_Core ${keyword} CGAL::CGAL )

message("libCGAL_Core is configured")

