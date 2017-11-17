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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_IO_PLY_READER_H
#define CGAL_IO_PLY_READER_H

#include <CGAL/IO/read_ply_points.h>

namespace CGAL{

  namespace internal
  {
    template <typename Integer, class Polygon_3, class Color_rgb>
    bool
    read_PLY_faces (std::istream& in,
                    PLY::PLY_reader& reader,
                    std::vector< Polygon_3 >& polygons,
                    std::vector< Color_rgb >& fcolors,
                    const char* vertex_indices_tag)
    {
      std::size_t faces_read = 0;

      bool has_colors = false;
      std::string rtag = "r", gtag = "g", btag = "b";
      if ((reader.does_tag_exist<boost::uint8_t>("red") || reader.does_tag_exist<boost::uint8_t>("r")) &&
          (reader.does_tag_exist<boost::uint8_t>("green") || reader.does_tag_exist<boost::uint8_t>("g")) &&
          (reader.does_tag_exist<boost::uint8_t>("blue") || reader.does_tag_exist<boost::uint8_t>("b")))
      {
        has_colors = true;
        if (reader.does_tag_exist<boost::uint8_t>("red"))
        {
          rtag = "red"; gtag = "green"; btag = "blue";
        }
      }
      
      while (!(in.eof()) && faces_read < reader.m_nb_faces)
      {
        for (std::size_t i = 0; i < reader.readers().size (); ++ i)
          reader.readers()[i]->get (in);

        cpp11::tuple<std::vector<Integer>, boost::uint8_t, boost::uint8_t, boost::uint8_t> new_face; 

        if (has_colors)
        {
          PLY::process_properties (reader, new_face,
                                   std::make_pair (CGAL::make_nth_of_tuple_property_map<0>(new_face),
                                                   PLY_property<std::vector<Integer> >(vertex_indices_tag)),
                                   std::make_pair (CGAL::make_nth_of_tuple_property_map<1>(new_face),
                                                   PLY_property<boost::uint8_t>(rtag.c_str())),
                                   std::make_pair (CGAL::make_nth_of_tuple_property_map<2>(new_face),
                                                   PLY_property<boost::uint8_t>(gtag.c_str())),
                                   std::make_pair (CGAL::make_nth_of_tuple_property_map<3>(new_face),
                                                   PLY_property<boost::uint8_t>(btag.c_str())));

          fcolors.push_back (Color_rgb (get<1>(new_face), get<2>(new_face), get<3>(new_face)));
        }
        else
          PLY::process_properties (reader, new_face,
                                   std::make_pair (CGAL::make_nth_of_tuple_property_map<0>(new_face),
                                                   PLY_property<std::vector<Integer> >(vertex_indices_tag)));

        polygons.push_back (Polygon_3(get<0>(new_face).size()));
        for (std::size_t i = 0; i < get<0>(new_face).size(); ++ i)
          polygons.back()[i] = std::size_t(get<0>(new_face)[i]);

        ++ faces_read;
      }

      return !in.bad();
    }

  }


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

    std::vector<CGAL::Color> dummy;

    reader.read_faces();

    if (reader.does_tag_exist<std::vector<boost::int32_t> > ("vertex_indices"))
      return internal::read_PLY_faces<boost::int32_t> (in, reader, polygons, dummy, "vertex_indices");

    if (reader.does_tag_exist<std::vector<boost::uint32_t> > ("vertex_indices"))
      return internal::read_PLY_faces<boost::uint32_t> (in, reader, polygons, dummy, "vertex_indices");

    if (reader.does_tag_exist<std::vector<boost::int32_t> > ("vertex_index"))
      return internal::read_PLY_faces<boost::int32_t> (in, reader, polygons, dummy, "vertex_index");

    if (reader.does_tag_exist<std::vector<boost::uint32_t> > ("vertex_index"))
      return internal::read_PLY_faces<boost::uint32_t> (in, reader, polygons, dummy, "vertex_index");

    std::cerr << "Error: can't find vertex indices in PLY input" << std::endl;
    return false;
  }

  template <class Point_3, class Polygon_3, class Color_rgb>
  bool
  read_PLY( std::istream& in,
            std::vector< Point_3 >& points,
            std::vector< Polygon_3 >& polygons,
            std::vector<Color_rgb>& fcolors,
            std::vector<Color_rgb>& vcolors,
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

    bool has_colors = false;
    std::string rtag = "r", gtag = "g", btag = "b";
    if ((reader.does_tag_exist<boost::uint8_t>("red") || reader.does_tag_exist<boost::uint8_t>("r")) &&
        (reader.does_tag_exist<boost::uint8_t>("green") || reader.does_tag_exist<boost::uint8_t>("g")) &&
        (reader.does_tag_exist<boost::uint8_t>("blue") || reader.does_tag_exist<boost::uint8_t>("b")))
    {
      has_colors = true;
      if (reader.does_tag_exist<boost::uint8_t>("red"))
      {
        rtag = "red"; gtag = "green"; btag = "blue";
      }
    }

    while (!(in.eof()) && points_read < reader.m_nb_points)
    {
      for (std::size_t i = 0; i < reader.readers().size (); ++ i)
        reader.readers()[i]->get (in);

      cpp11::tuple<Point_3, boost::uint8_t, boost::uint8_t, boost::uint8_t> new_vertex;

      if (has_colors)
      {
        internal::PLY::process_properties (reader, new_vertex,
                                           make_ply_point_reader (CGAL::make_nth_of_tuple_property_map<0>(new_vertex)),
                                           std::make_pair (CGAL::make_nth_of_tuple_property_map<1>(new_vertex),
                                                           PLY_property<boost::uint8_t>(rtag.c_str())),
                                           std::make_pair (CGAL::make_nth_of_tuple_property_map<2>(new_vertex),
                                                           PLY_property<boost::uint8_t>(gtag.c_str())),
                                           std::make_pair (CGAL::make_nth_of_tuple_property_map<3>(new_vertex),
                                                           PLY_property<boost::uint8_t>(btag.c_str())));

        vcolors.push_back (Color_rgb (get<1>(new_vertex), get<2>(new_vertex), get<3>(new_vertex)));
      }
      else
        internal::PLY::process_properties (reader, new_vertex,
                                           make_ply_point_reader (CGAL::make_nth_of_tuple_property_map<0>(new_vertex)));
      
      points.push_back (get<0>(new_vertex));
      
      ++ points_read;
    }

    if (points_read != reader.m_nb_points)
      return false;

    reader.read_faces();
    
    if (reader.does_tag_exist<std::vector<boost::int32_t> > ("vertex_indices"))
      return internal::read_PLY_faces<boost::int32_t> (in, reader, polygons, fcolors, "vertex_indices");

    if (reader.does_tag_exist<std::vector<boost::uint32_t> > ("vertex_indices"))
      return internal::read_PLY_faces<boost::uint32_t> (in, reader, polygons, fcolors, "vertex_indices");

    if (reader.does_tag_exist<std::vector<boost::int32_t> > ("vertex_index"))
      return internal::read_PLY_faces<boost::int32_t> (in, reader, polygons, fcolors, "vertex_index");

    if (reader.does_tag_exist<std::vector<boost::uint32_t> > ("vertex_index"))
      return internal::read_PLY_faces<boost::uint32_t> (in, reader, polygons, fcolors, "vertex_index");

    std::cerr << "Error: can't find vertex indices in PLY input" << std::endl;
    return false;
  }


} // namespace CGAL

#endif // CGAL_IO_PLY_READER_H
