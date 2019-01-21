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
                    internal::PLY::PLY_element& element,
                    std::vector< Polygon_3 >& polygons,
                    std::vector< Color_rgb >& fcolors,
                    const char* vertex_indices_tag)
    {
      bool has_colors = false;
      std::string rtag = "r", gtag = "g", btag = "b";
      if ((element.has_property<boost::uint8_t>("red") || element.has_property<boost::uint8_t>("r")) &&
          (element.has_property<boost::uint8_t>("green") || element.has_property<boost::uint8_t>("g")) &&
          (element.has_property<boost::uint8_t>("blue") || element.has_property<boost::uint8_t>("b")))
      {
        has_colors = true;
        if (element.has_property<boost::uint8_t>("red"))
        {
          rtag = "red"; gtag = "green"; btag = "blue";
        }
      }
      
      for (std::size_t j = 0; j < element.number_of_items(); ++ j)
      {
        for (std::size_t k = 0; k < element.number_of_properties(); ++ k)
        {
          internal::PLY::PLY_read_number* property = element.property(k);
          property->get (in);

          if (in.fail())
            return false;
        }

        cpp11::tuple<std::vector<Integer>, boost::uint8_t, boost::uint8_t, boost::uint8_t> new_face; 

        if (has_colors)
        {
          PLY::process_properties (element, new_face,
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
          PLY::process_properties (element, new_face,
                                   std::make_pair (CGAL::make_nth_of_tuple_property_map<0>(new_face),
                                                   PLY_property<std::vector<Integer> >(vertex_indices_tag)));

        polygons.push_back (Polygon_3(get<0>(new_face).size()));
        for (std::size_t i = 0; i < get<0>(new_face).size(); ++ i)
          polygons.back()[i] = std::size_t(get<0>(new_face)[i]);
      }

      return true;
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
    {
      in.setstate(std::ios::failbit);
      return false;
    }
  
    for (std::size_t i = 0; i < reader.number_of_elements(); ++ i)
    {
      internal::PLY::PLY_element& element = reader.element(i);

      if (element.name() == "vertex" || element.name() == "vertices")
      {
        for (std::size_t j = 0; j < element.number_of_items(); ++ j)
        {
          for (std::size_t k = 0; k < element.number_of_properties(); ++ k)
          {
            internal::PLY::PLY_read_number* property = element.property(k);
            property->get (in);

            if (in.fail())
              return false;
          }

          Point_3 new_vertex;

          internal::PLY::process_properties (element, new_vertex,
                                             make_ply_point_reader (CGAL::Identity_property_map<Point_3>()));
      
          points.push_back (get<0>(new_vertex));
        }
      }
      else if (element.name() == "face" || element.name() == "faces")
      {
        std::vector<CGAL::Color> dummy;

        if (element.has_property<std::vector<boost::int32_t> > ("vertex_indices"))
          internal::read_PLY_faces<boost::int32_t> (in, element, polygons, dummy, "vertex_indices");
        else if (element.has_property<std::vector<boost::uint32_t> > ("vertex_indices"))
          internal::read_PLY_faces<boost::uint32_t> (in, element, polygons, dummy, "vertex_indices");
        else if (element.has_property<std::vector<boost::int32_t> > ("vertex_index"))
          internal::read_PLY_faces<boost::int32_t> (in, element, polygons, dummy, "vertex_index");
        else if (element.has_property<std::vector<boost::uint32_t> > ("vertex_index"))
          internal::read_PLY_faces<boost::uint32_t> (in, element, polygons, dummy, "vertex_index");
        else
        {
          std::cerr << "Error: can't find vertex indices in PLY input" << std::endl;
          return false;
        }
      }
      else // Read other elements and ignore
      {
        for (std::size_t j = 0; j < element.number_of_items(); ++ j)
        {
          for (std::size_t k = 0; k < element.number_of_properties(); ++ k)
          {
            internal::PLY::PLY_read_number* property = element.property(k);
            property->get (in);

            if (in.fail())
              return false;
          }
        }
      }
    }

    return true;
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
    {
      in.setstate(std::ios::failbit);
      return false;
    }

    for (std::size_t i = 0; i < reader.number_of_elements(); ++ i)
    {
      internal::PLY::PLY_element& element = reader.element(i);

      if (element.name() == "vertex" || element.name() == "vertices")
      {
        bool has_colors = false;
        std::string rtag = "r", gtag = "g", btag = "b";
        if ((element.has_property<boost::uint8_t>("red") || element.has_property<boost::uint8_t>("r")) &&
            (element.has_property<boost::uint8_t>("green") || element.has_property<boost::uint8_t>("g")) &&
            (element.has_property<boost::uint8_t>("blue") || element.has_property<boost::uint8_t>("b")))
        {
          has_colors = true;
          if (element.has_property<boost::uint8_t>("red"))
          {
            rtag = "red"; gtag = "green"; btag = "blue";
          }
        }

        for (std::size_t j = 0; j < element.number_of_items(); ++ j)
        {
          for (std::size_t k = 0; k < element.number_of_properties(); ++ k)
          {
            internal::PLY::PLY_read_number* property = element.property(k);
            property->get (in);

            if (in.fail())
              return false;
          }

          cpp11::tuple<Point_3, boost::uint8_t, boost::uint8_t, boost::uint8_t> new_vertex;

          if (has_colors)
          {
            internal::PLY::process_properties (element, new_vertex,
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
            internal::PLY::process_properties (element, new_vertex,
                                               make_ply_point_reader (CGAL::make_nth_of_tuple_property_map<0>(new_vertex)));

          points.push_back (get<0>(new_vertex));
        }
      }
      else if (element.name() == "face" || element.name() == "faces")
      {
        if (element.has_property<std::vector<boost::int32_t> > ("vertex_indices"))
          internal::read_PLY_faces<boost::int32_t> (in, element, polygons, fcolors, "vertex_indices");
        else if (element.has_property<std::vector<boost::uint32_t> > ("vertex_indices"))
          internal::read_PLY_faces<boost::uint32_t> (in, element, polygons, fcolors, "vertex_indices");
        else if (element.has_property<std::vector<boost::int32_t> > ("vertex_index"))
          internal::read_PLY_faces<boost::int32_t> (in, element, polygons, fcolors, "vertex_index");
        else if (element.has_property<std::vector<boost::uint32_t> > ("vertex_index"))
          internal::read_PLY_faces<boost::uint32_t> (in, element, polygons, fcolors, "vertex_index");
        else
        {
          std::cerr << "Error: can't find vertex indices in PLY input" << std::endl;
          return false;
        }
      }
      else // Read other elements and ignore
      {
        for (std::size_t j = 0; j < element.number_of_items(); ++ j)
        {
          for (std::size_t k = 0; k < element.number_of_properties(); ++ k)
          {
            internal::PLY::PLY_read_number* property = element.property(k);
            property->get (in);

            if (in.fail())
              return false;
          }
        }
      }
    }

    return true;
  }


} // namespace CGAL

#endif // CGAL_IO_PLY_READER_H
