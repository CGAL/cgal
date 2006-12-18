// Copyright (c) 2006 Geometry Factory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s): Andreas Fabri <andreas.fabri@geometryfactory.com>, Fernando Cacciola <fernando.cacciola@gmail.com>

#ifndef CGAL_BOOST_GRAPH_BGL_PROPERTIES_H
#define CGAL_BOOST_GRAPH_BGL_PROPERTIES_H

#include <boost/property_map.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

enum vertex_is_border_t      { vertex_is_border } ;
enum vertex_point_t          { vertex_point     } ;
enum vertex_external_index_t { vertex_external_index } ;
enum edge_is_border_t        { edge_is_border   } ;
enum edge_external_index_t   { edge_external_index } ;

CGAL_END_NAMESPACE

#endif // CGAL_BOOST_GRAPH_BGL_PROPERTIES_H
