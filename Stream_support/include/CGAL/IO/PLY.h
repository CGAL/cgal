// Copyright (c) 2015  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_IO_PLY_H
#define CGAL_IO_PLY_H

#include <CGAL/IO/PLY/PLY_reader.h>
#include <CGAL/IO/PLY/PLY_writer.h>
#include <CGAL/property_map.h>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

template <class Point_3, class Polygon_3, class Color_rgb>
bool read_PLY(std::istream& in,
              std::vector<Point_3>& points,
              std::vector<Polygon_3>& polygons,
              std::vector<std::pair<unsigned int, unsigned int> >& hedges,
              std::vector<Color_rgb>& fcolors,
              std::vector<Color_rgb>& vcolors,
              std::vector<std::pair<float, float> >& huvs,
              bool /* verbose */ = false)
{
  if(!in)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }
  internal::PLY::PLY_reader reader;

  if(!(reader.init(in)))
  {
    in.setstate(std::ios::failbit);
    return false;
  }
  for(std::size_t i = 0; i < reader.number_of_elements(); ++ i)
  {
    internal::PLY::PLY_element& element = reader.element(i);

    if(element.name() == "vertex" || element.name() == "vertices")
    {
      bool has_colors = false;
      std::string rtag = "r", gtag = "g", btag = "b";
      if((element.has_property<boost::uint8_t>("red") || element.has_property<boost::uint8_t>("r")) &&
         (element.has_property<boost::uint8_t>("green") || element.has_property<boost::uint8_t>("g")) &&
         (element.has_property<boost::uint8_t>("blue") || element.has_property<boost::uint8_t>("b")))
      {
        has_colors = true;
        if(element.has_property<boost::uint8_t>("red"))
        {
          rtag = "red"; gtag = "green"; btag = "blue";
        }
      }

      for(std::size_t j = 0; j < element.number_of_items(); ++ j)
      {
        for(std::size_t k = 0; k < element.number_of_properties(); ++ k)
        {
          internal::PLY::PLY_read_number* property = element.property(k);
          property->get(in);

          if(in.fail())
            return false;
        }

        std::tuple<Point_3, boost::uint8_t, boost::uint8_t, boost::uint8_t> new_vertex;

        if(has_colors)
        {
          internal::PLY::process_properties(element, new_vertex,
                                            make_ply_point_reader(CGAL::make_nth_of_tuple_property_map<0>(new_vertex)),
                                            std::make_pair(CGAL::make_nth_of_tuple_property_map<1>(new_vertex),
                                                           PLY_property<boost::uint8_t>(rtag.c_str())),
                                            std::make_pair(CGAL::make_nth_of_tuple_property_map<2>(new_vertex),
                                                           PLY_property<boost::uint8_t>(gtag.c_str())),
                                            std::make_pair(CGAL::make_nth_of_tuple_property_map<3>(new_vertex),
                                                           PLY_property<boost::uint8_t>(btag.c_str())));

          vcolors.push_back(Color_rgb(get<1>(new_vertex), get<2>(new_vertex), get<3>(new_vertex)));
        }
        else
          internal::PLY::process_properties(element, new_vertex,
                                            make_ply_point_reader(CGAL::make_nth_of_tuple_property_map<0>(new_vertex)));

        points.push_back(get<0>(new_vertex));
      }
    }
    else if(element.name() == "face" || element.name() == "faces")
    {
      if(element.has_property<std::vector<boost::int32_t> >("vertex_indices"))
        internal::read_PLY_faces<boost::int32_t>(in, element, polygons, fcolors, "vertex_indices");
      else if(element.has_property<std::vector<boost::uint32_t> >("vertex_indices"))
        internal::read_PLY_faces<boost::uint32_t>(in, element, polygons, fcolors, "vertex_indices");
      else if(element.has_property<std::vector<boost::int32_t> >("vertex_index"))
        internal::read_PLY_faces<boost::int32_t>(in, element, polygons, fcolors, "vertex_index");
      else if(element.has_property<std::vector<boost::uint32_t> >("vertex_index"))
        internal::read_PLY_faces<boost::uint32_t>(in, element, polygons, fcolors, "vertex_index");
      else
      {
        std::cerr << "Error: can't find vertex indices in PLY input" << std::endl;
        return false;
      }
    }
    else if(element.name() == "halfedge" )
    {
      bool has_uv = false;
      std::string stag = "source", ttag = "target", utag = "u", vtag = "v";
      if( element.has_property<unsigned int>("source") &&
          element.has_property<unsigned int>("target") &&
          element.has_property<float>("u") &&
          element.has_property<float>("v"))
      {
        has_uv = true;
      }
      cpp11::tuple<unsigned int, unsigned int, float, float, float>  new_hedge;
      for(std::size_t j = 0; j < element.number_of_items(); ++ j)
      {
        for(std::size_t k = 0; k < element.number_of_properties(); ++ k)
        {
          internal::PLY::PLY_read_number* property = element.property(k);
          property->get(in);

          if(in.eof())
            return false;
        }

        if(has_uv)
        {
          internal::PLY::process_properties(element, new_hedge,
                                            std::make_pair(CGAL::make_nth_of_tuple_property_map<0>(new_hedge),
                                                           PLY_property<unsigned int>(stag.c_str())),
                                            std::make_pair(CGAL::make_nth_of_tuple_property_map<1>(new_hedge),
                                                           PLY_property<unsigned int>(ttag.c_str())),
                                            std::make_pair(CGAL::make_nth_of_tuple_property_map<2>(new_hedge),
                                                           PLY_property<float>(utag.c_str())),
                                            std::make_pair(CGAL::make_nth_of_tuple_property_map<3>(new_hedge),
                                                           PLY_property<float>(vtag.c_str())));
          hedges.push_back(std::make_pair(get<0>(new_hedge), get<1>(new_hedge)));
          huvs.push_back(std::make_pair(get<2>(new_hedge), get<3>(new_hedge)));
        }
        else
        {
          internal::PLY::process_properties(element, new_hedge,
                                            std::make_pair(CGAL::make_nth_of_tuple_property_map<0>(new_hedge),
                                                           PLY_property<unsigned int>(stag.c_str())),
                                            std::make_pair(CGAL::make_nth_of_tuple_property_map<1>(new_hedge),
                                                           PLY_property<unsigned int>(ttag.c_str()))
                                            );
        }
      }
    }
    else // Read other elements and ignore
    {
      for(std::size_t j = 0; j < element.number_of_items(); ++ j)
      {
        for(std::size_t k = 0; k < element.number_of_properties(); ++ k)
        {
          internal::PLY::PLY_read_number* property = element.property(k);
          property->get(in);
          if(in.fail())
            return false;
        }
      }
    }
  }
  return !in.bad();
}

template <class Point_3, class Polygon_3, class Color_rgb>
bool read_PLY(std::istream& in,
              std::vector< Point_3 >& points,
              std::vector< Polygon_3 >& polygons,
              std::vector<Color_rgb>& fcolors,
              std::vector<Color_rgb>& vcolors,
              bool /* verbose */ = false)
{
  std::vector<std::pair<unsigned int, unsigned int> > dummy_pui;
  std::vector<std::pair<float, float> > dummy_pf;

  return read_PLY<Point_3, Polygon_3, Color_rgb>(in, points, polygons, dummy_pui, fcolors, vcolors, dummy_pf);
}

/*!
 * \ingroup IOstreamFunctions
 *
 * reads the content of `in` into `points` and `polygons`, in the PLY format.
 *
 * \see \ref IOStreamPLY
 */
template <class Point_3, class Polygon_3>
bool
read_PLY(std::istream& in,
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

  if(!(reader.init(in)))
  {
    in.setstate(std::ios::failbit);
    return false;
  }

  for(std::size_t i = 0; i < reader.number_of_elements(); ++ i)
  {
    internal::PLY::PLY_element& element = reader.element(i);

    if(element.name() == "vertex" || element.name() == "vertices")
    {
      for(std::size_t j = 0; j < element.number_of_items(); ++ j)
      {
        for(std::size_t k = 0; k < element.number_of_properties(); ++ k)
        {
          internal::PLY::PLY_read_number* property = element.property(k);
          property->get(in);

          if(in.fail())
            return false;
        }

        Point_3 new_vertex;

        internal::PLY::process_properties(element, new_vertex,
                                          make_ply_point_reader(CGAL::Identity_property_map<Point_3>()));

        points.push_back(get<0>(new_vertex));
      }
    }
    else if(element.name() == "face" || element.name() == "faces")
    {
      std::vector<CGAL::Color> dummy;

      if(element.has_property<std::vector<boost::int32_t> >("vertex_indices"))
        internal::read_PLY_faces<boost::int32_t>(in, element, polygons, dummy, "vertex_indices");
      else if(element.has_property<std::vector<boost::uint32_t> >("vertex_indices"))
        internal::read_PLY_faces<boost::uint32_t>(in, element, polygons, dummy, "vertex_indices");
      else if(element.has_property<std::vector<boost::int32_t> >("vertex_index"))
        internal::read_PLY_faces<boost::int32_t>(in, element, polygons, dummy, "vertex_index");
      else if(element.has_property<std::vector<boost::uint32_t> >("vertex_index"))
        internal::read_PLY_faces<boost::uint32_t>(in, element, polygons, dummy, "vertex_index");
      else
      {
        std::cerr << "Error: can't find vertex indices in PLY input" << std::endl;
        return false;
      }
    }
    else // Read other elements and ignore
    {
      for(std::size_t j = 0; j < element.number_of_items(); ++ j)
      {
        for(std::size_t k = 0; k < element.number_of_properties(); ++ k)
        {
          internal::PLY::PLY_read_number* property = element.property(k);
          property->get(in);

          if(in.fail())
            return false;
        }
      }
    }
  }

  return !in.bad();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/*!
 * \ingroup IOstreamFunctions
 *
 * writes the content of `points` and `polygons` in `out`, in the OFF format.
 *
 * \see \ref IOStreamOFF
 */
template <class Point_3, class Polygon_3>
bool write_PLY(std::ostream& out,
               const std::vector< Point_3 >& points,
               const std::vector< Polygon_3 >& polygons,
               bool /* verbose */ = false)
{
  if(!out)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // Write header
  out << "ply" << std::endl
      << ((get_mode(out) == IO::BINARY) ? "format binary_little_endian 1.0" : "format ascii 1.0") << std::endl
      << "comment Generated by the CGAL library" << std::endl
      << "element vertex " << points.size() << std::endl;

  internal::PLY::output_property_header(out, make_ply_point_writer (CGAL::Identity_property_map<Point_3>()));

  out << "element face " << polygons.size() << std::endl;

  internal::PLY::output_property_header(out, std::make_pair(CGAL::Identity_property_map<Polygon_3>(),
                                                            PLY_property<std::vector<int> >("vertex_indices")));

  out << "end_header" << std::endl;

  for (std::size_t i = 0; i < points.size(); ++ i)
    internal::PLY::output_properties(out, points.begin() + i,
                                     make_ply_point_writer (CGAL::Identity_property_map<Point_3>()));

  for (std::size_t i = 0; i < polygons.size(); ++ i)
    internal::PLY::output_properties(out, polygons.begin() + i,
                                     std::make_pair(CGAL::Identity_property_map<Polygon_3>(),
                                                    PLY_property<std::vector<int> >("vertex_indices")));

  return out.good();
}

template <class SurfaceMesh>
bool write_PLY(std::ostream& out,
               const SurfaceMesh& mesh,
               bool /* verbose */ = false)
{
  typedef typename boost::graph_traits<SurfaceMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<SurfaceMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<SurfaceMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_map<SurfaceMesh, boost::vertex_point_t>::type::value_type Point_3;
  typedef typename SurfaceMesh::template Property_map<halfedge_descriptor,std::pair<float, float> > UV_map;

  UV_map h_uv;
  bool has_texture;
  boost::tie(h_uv, has_texture) = mesh.template property_map<halfedge_descriptor,std::pair<float, float> >("h:uv");

  if(!out)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // Write header
  out << "ply" << std::endl
      << ((get_mode(out) == IO::BINARY) ? "format binary_little_endian 1.0" : "format ascii 1.0") << std::endl
      << "comment Generated by the CGAL library" << std::endl
      << "element vertex " << num_vertices(mesh) << std::endl;

  internal::PLY::output_property_header(out,
                                        make_ply_point_writer (CGAL::Identity_property_map<Point_3>()));

  out << "element face " << num_faces(mesh) << std::endl;

  internal::PLY::output_property_header(out,
                                        std::make_pair(CGAL::Identity_property_map<std::vector<std::size_t> >(),
                                                       PLY_property<std::vector<int> >("vertex_indices")));

  if(has_texture)
  {
    out << "element halfedge " << num_halfedges(mesh) << std::endl;

    internal::PLY::output_property_header(out,
                                          std::make_pair(CGAL::Identity_property_map<std::size_t >(),
                                                         PLY_property<unsigned int >("source")));

    internal::PLY::output_property_header(out,
                                          std::make_pair(CGAL::Identity_property_map<std::size_t >(),
                                                         PLY_property<unsigned int >("target")));
    internal::PLY::output_property_header(out,
                                          std::make_tuple (h_uv,
                                                           PLY_property<float>("u"),
                                                           PLY_property<float>("v")));
  }
  out << "end_header" << std::endl;

  for(vertex_descriptor vd : vertices(mesh))
  {
    Point_3 p = get(get(CGAL::vertex_point, mesh), vd);
    internal::PLY::output_properties(out, &p,
                                     make_ply_point_writer (CGAL::Identity_property_map<Point_3>()));
  }

  std::vector<std::size_t> polygon;
  for(face_descriptor fd : faces(mesh))
  {
    polygon.clear();
    for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, mesh), mesh))
      polygon.push_back(get(get(boost::vertex_index, mesh), target(hd,mesh)));

    internal::PLY::output_properties(out, &polygon,
                                     std::make_pair(CGAL::Identity_property_map<std::vector<std::size_t> >(),
                                                    PLY_property<std::vector<int> >("vertex_indices")));
  }

  if(has_texture)
  {
    for(halfedge_descriptor hd : halfedges(mesh))
    {
      typedef std::tuple<unsigned int, unsigned int, float, float> Super_tuple;
      Super_tuple t = std::make_tuple(source(hd, mesh),target(hd, mesh),
                                      h_uv[hd].first,
                                      h_uv[hd].second);

      internal::PLY::output_properties(out, &t,
                                       std::make_pair(Nth_of_tuple_property_map<0,Super_tuple>(),
                                                      PLY_property<unsigned int >("source")),
                                       std::make_pair(Nth_of_tuple_property_map<1,Super_tuple>(),
                                                      PLY_property<unsigned int >("target")),
                                       std::make_pair(Nth_of_tuple_property_map<2,Super_tuple>(),
                                                      PLY_property<float>("u")),
                                       std::make_pair(Nth_of_tuple_property_map<3,Super_tuple>(),
                                                      PLY_property<float>("v")));
    }
  }

  return out.good();
}

} // namespace CGAL

#endif // CGAL_IO_PLY_H
