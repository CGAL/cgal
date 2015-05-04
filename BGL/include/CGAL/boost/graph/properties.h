// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
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
//
// Author(s)     : Andreas Fabri, Fernando Cacciola


#ifndef CGAL_BOOST_GRAPH_BGL_PROPERTIES_H
#define CGAL_BOOST_GRAPH_BGL_PROPERTIES_H

#include <CGAL/property_map.h>
#include <boost/graph/properties.hpp>

#include <CGAL/basic.h>
#include <string>

/// Boost Namespace
namespace boost {

/// \ingroup PkgBGLProperties
/// @{

/// A property tag which refers to the geometric embedding property
/// of a vertex of a \ref HalfedgeGraph.
enum vertex_point_t          { vertex_point          };
enum vertex_external_index_t { vertex_external_index } ;

/// A property tag which refers to the property
/// of a halfedge of being a border halfedge.
enum edge_external_index_t   { edge_external_index   } ;

/// A property tag which identifies the *index* property of
/// a halfedge of a \ref HalfedgeGraph.
enum halfedge_index_t        { halfedge_index        };
enum halfedge_external_index_t   { halfedge_external_index   } ;

/// A property tag which identifies the *index* property of
/// a face of a \ref FaceGraph.
enum face_index_t            { face_index            };
enum face_external_index_t   { face_external_index   } ;

  template <typename T>
  struct vertex_property_t
  {
    vertex_property_t(const std::string s, const T& t = T())
      : s(s), t(t)
    {}
    std::string s;
    T t;
        
  };

/// @}

// Introduce those two tags so we can use BOOST_INSTALL_PROPERTY
// macro. This is dangerous because we now rely on implementation
// details.
struct halfedge_property_tag { };
struct face_property_tag { };

BOOST_INSTALL_PROPERTY(vertex, point);
BOOST_INSTALL_PROPERTY(vertex, external_index);
BOOST_INSTALL_PROPERTY(halfedge, external_index);
BOOST_INSTALL_PROPERTY(edge, external_index);
BOOST_INSTALL_PROPERTY(face, index);
BOOST_INSTALL_PROPERTY(face, external_index);
} // boost

namespace CGAL {
using boost::vertex_point_t;
using boost::vertex_point;
using boost::vertex_external_index_t;
using boost::vertex_external_index;
using boost::halfedge_index_t;
using boost::halfedge_index;
using boost::halfedge_external_index_t;
using boost::halfedge_external_index;
using boost::edge_external_index_t;
using boost::edge_external_index;
using boost::face_index_t;
using boost::face_index;
using boost::face_external_index_t;
using boost::face_external_index;
using boost::vertex_property_t;
} // CGAL

#endif // CGAL_BOOST_GRAPH_BGL_PROPERTIES_H
