set(list_of_whitelisted_headers_txt [=[
  CGAL/Linear_cell_complex_constructors.h
  CGAL/CGAL_Ipelet_base.h
]=])

separate_arguments(list_of_whitelisted_headers UNIX_COMMAND ${list_of_whitelisted_headers_txt})
