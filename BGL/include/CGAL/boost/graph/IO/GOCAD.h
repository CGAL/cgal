// Copyright (c) 2015  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//                 Mael Rouxel-Labbé

#ifndef CGAL_BGL_IO_GOCAD_H
#define CGAL_BGL_IO_GOCAD_H

#include <CGAL/boost/graph/IO/Generic_facegraph_builder.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/container/flat_map.hpp>

#include <fstream>
#include <iostream>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

namespace IO {
namespace internal {

// Use CRTP to gain access to the protected members without getters/setters.
template <typename FaceGraph, typename Point>
class GOCAD_builder
  : public Generic_facegraph_builder<FaceGraph, Point, GOCAD_builder<FaceGraph, Point> >
{
  typedef GOCAD_builder<FaceGraph, Point>                                       Self;
  typedef Generic_facegraph_builder<FaceGraph, Point, Self>                     Base;

  typedef typename Base::Point_container                                        Point_container;
  typedef typename Base::Face                                                   Face;
  typedef typename Base::Face_container                                         Face_container;

public:
  GOCAD_builder(std::istream& is_) : Base(is_) { }

  // Implementation of the two functions required by the generic builder
  bool read(std::istream& input,
            Point_container& points,
            Face_container& faces)
  {
    int offset = 0;
    char c;
    std::string s, tface("TFACE");
    int i,j,k;
    Point p;
    bool vertices_read = false;

    while(input >> s)
    {
      if(s == tface)
        break;

      std::string::size_type idx;
      if((idx = s.find("name")) != std::string::npos)
      {
        std::istringstream str(s.substr(idx + 5));
        str >> this->name;
      }

      if((idx = s.find("color")) != std::string::npos)
      {
        std::istringstream str(s.substr(idx + 6));
        str >> this->color;
      }
    }
    std::getline(input, s);

    while(input.get(c))
    {
      if((c == 'V') || (c == 'P'))
      {
        input >> s >> i >> p; // @fixme check for failure
        if(!vertices_read)
        {
          vertices_read = true;
          offset -= i; // Some files start with index 0 others with 1
        }

        points.push_back(p);
      }
      else if(vertices_read && (c == 'T'))
      {
        input >> c >> c >> c >>  i >> j >> k;
        typename Base::Face new_face(3);
        new_face[0] = offset+i;
        new_face[1] = offset+j;
        new_face[2] = offset+k;
        faces.push_back(new_face);
      }
      else if(c == 'E')
      {
        break;
      }

      std::getline(input, s);
    }

    return true;
  }
};

} // namespace internal
} // namespace IO

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the TS format.
  `name` and `color` will be filled according to the values contained in the file.

  \pre The data must represent a 2-manifold

  \attention The graph `g` is not cleared, and the data from the stream are added.

  \see \ref IOStreamGocad
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_GOCAD(std::istream& in,
                std::string& name,
                std::string& color,
                FaceGraph& g,
                const CGAL_BGL_NP_CLASS& np)
{
  typedef typename CGAL::GetVertexPointMap<FaceGraph, CGAL_BGL_NP_CLASS>::type VPM;
  typedef typename boost::property_traits<VPM>::value_type                     Point;

  IO::internal::GOCAD_builder<FaceGraph, Point> builder(in);
  if(!builder(g, np))
    return false;

  name = builder.name;
  color = builder.color;

  return g.is_valid(); // @fixme keep validity check?
}

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the TS format.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \pre The data must represent a 2-manifold

  \see \ref IOStreamGOCAD
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_GOCAD(const char* fname, FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  std::ifstream in(fname);
  std::string unused_name;
  std::string unused_color;

  return read_GOCAD(in, unused_name, unused_color, g, np);
}

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_GOCAD(const std::string& fname, FaceGraph& g, CGAL_BGL_NP_CLASS np)
{
  return read_GOCAD(fname.c_str(), g, np);
}

template <typename FaceGraph>
bool read_GOCAD(std::istream& is, FaceGraph& g) { return read_GOCAD(is, g, parameters::all_default()); }
template <typename FaceGraph>
bool read_GOCAD(const char* fname, FaceGraph& g) { return read_GOCAD(fname, g, parameters::all_default()); }
template <typename FaceGraph>
bool read_GOCAD(const std::string& fname, FaceGraph& g) { return read_GOCAD(fname, g, parameters::all_default()); }

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the TS format into `os`. `fname` is the
  mandatory name that will be assigned to `g`in the file.

  \see \ref IOStreamGocad
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(std::ostream& os,
                 const char* fname,
                 const FaceGraph& g,
                 const CGAL_BGL_NP_CLASS& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor      face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type   vertices_size_type;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typename CGAL::GetVertexPointMap<FaceGraph, CGAL_BGL_NP_CLASS>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, g));

  if(!os.good())
    return false;

  os << "GOCAD TSurf 1\n"
        "HEADER {\n"
        "name:";
  os << fname << std::endl;
  os << "*border:on\n"
        "*border*bstone:on\n"
        "}\n"
        "GOCAD_ORIGINAL_COORDINATE_SYSTEM\n"
        "NAME Default\n"
        "AXIS_NAME \"X\" \"Y\" \"Z\"\n"
        "AXIS_UNIT \"m\" \"m\" \"m\"\n"
        "ZPOSITIVE Elevation\n"
        "END_ORIGINAL_COORDINATE_SYSTEM\n"
        "TFACE\n";

  os.precision(16);

  boost::container::flat_map<vertex_descriptor, vertices_size_type> reindex;

  vertices_size_type i = 0;
  for(const vertex_descriptor v : vertices(g))
  {
    os << "VRTX " << i << " " << get(vpm, v) << "\n";
    reindex[v] = i++;
  }

  for(const face_descriptor f : faces(g))
  {
    halfedge_descriptor h = halfedge(f, g);
    os << "TRGL " << reindex[target(prev(h, g), g)] << " "
                  << reindex[target(h, g)] << " "
                  << reindex[target(next(h, g), g)] << "\n";
  }

  os << "END" << std::endl;

  return os.good();
}

/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the TS format into a file named `fname`.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \see \ref IOStreamGOCAD
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(const char* fname,
                 const FaceGraph& g,
                 const CGAL_BGL_NP_CLASS& np)
{
  std::ofstream out(fname);
  return write_GOCAD(out, fname, g, np);
}

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_GOCAD(const std::string& fname, const FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  return write_GOCAD(fname.c_str(), g, np);
}

template <typename FaceGraph>
bool write_GOCAD(const char* fname, const FaceGraph& g) { return write_GOCAD(fname, g, parameters::all_default()); }
template <typename FaceGraph>
bool write_GOCAD(const std::string& fname, const FaceGraph& g) { return write_GOCAD(fname, g, parameters::all_default()); }

} // namespace CGAL

#endif // CGAL_BGL_IO_GOCAD_H
