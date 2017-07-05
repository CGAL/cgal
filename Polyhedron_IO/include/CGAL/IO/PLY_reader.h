// Copyright (c) 2017 GeometryFactory
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_IO_PLY_READER_H
#define CGAL_IO_PLY_READER_H

#include <CGAL/license/Polyhedron.h>

#include <CGAL/IO/read_ply_points.h>

namespace CGAL{

  template <class Point_3, class Polygon_3>
  bool
  read_PLY( std::istream& in,
            std::vector< Point_3 >& points,
            std::vector< Polygon_3 >& polygons,
            bool /* verbose */ = false)
  {
    if(!in)
    {
      std::cerr << "Error: cannot open file" << std::endl;
      return false;
    }

    internal::PLY::PLY_reader reader;
  
    if (!(reader.init (in)))
      return false;
  
    std::size_t points_read = 0;
  
    while (!(in.eof()) && points_read < reader.m_nb_points)
    {
      for (std::size_t i = 0; i < reader.readers().size (); ++ i)
        reader.readers()[i]->get (in);

      Point_3 new_vertex;

      internal::PLY::process_properties (reader, new_vertex,
                                         make_ply_point_reader (CGAL::Identity_property_map<Point_3>()));
                                         

      points.push_back (new_vertex);
      
      ++ points_read;
    }

    if (points_read != reader.m_nb_points)
      return false;

    reader.read_faces();

    std::size_t faces_read = 0;
  
    while (!(in.eof()) && faces_read < reader.m_nb_faces)
    {

      for (std::size_t i = 0; i < reader.readers().size (); ++ i)
        reader.readers()[i]->get (in);

      std::vector<int> new_face; 

      internal::PLY::process_properties (reader, new_face,
                                         std::make_pair (CGAL::Identity_property_map<std::vector<int> >(),
                                                         PLY_property<std::vector<int> >("vertex_indices")));

      polygons.push_back (Polygon_3(new_face.size()));
      for (std::size_t i = 0; i < new_face.size(); ++ i)
        polygons.back()[i] = std::size_t(new_face[i]);

      ++ faces_read;
    }


    return in.good();
  }


} // namespace CGAL

#endif // CGAL_IO_PLY_READER_H
