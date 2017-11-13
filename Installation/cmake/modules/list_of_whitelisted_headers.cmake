set(list_of_whitelisted_headers_txt [=[
  CGAL/Linear_cell_complex_constructors.h
  CGAL/CGAL_Ipelet_base.h
  CGAL/IO/read_las_points.h
  CGAL/IO/write_las_points.h
  CGAL/IO/read_ply_points.h
  CGAL/IO/write_ply_points.h
]=])

separate_arguments(list_of_whitelisted_headers UNIX_COMMAND ${list_of_whitelisted_headers_txt})
