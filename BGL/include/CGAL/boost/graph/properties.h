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

namespace CGAL {

enum vertex_is_border_t      { vertex_is_border } ;
enum vertex_point_t          { vertex_point     } ;
enum vertex_external_index_t { vertex_external_index } ;
enum edge_is_border_t        { edge_is_border   } ;
enum edge_external_index_t   { edge_external_index } ;

} //namespace CGAL

#endif // CGAL_BOOST_GRAPH_BGL_PROPERTIES_H
