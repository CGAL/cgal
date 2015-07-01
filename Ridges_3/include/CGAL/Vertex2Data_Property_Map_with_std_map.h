// Copyright (c) 2015  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Marc Pouget and Frédéric Cazals

#ifndef CGAL_VERTEX2DATA_PROPERTY_MAP_WITH_STD_MAP_H
#define CGAL_VERTEX2DATA_PROPERTY_MAP_WITH_STD_MAP_H 

//---------------------------------------------------------------------------
//Vertex2Data_Property_Map_with_std_map
// defines models for Vertex2FTPropertyMap and Vertex2VectorPropertyMap
//--------------------------------------------------------------------------
template < class TriangulatedSurfaceMesh >
class Vertex2Data_Property_Map_with_std_map 
{
 public:
  typedef typename TriangulatedSurfaceMesh::Traits::FT        FT;
  typedef typename TriangulatedSurfaceMesh::Traits::Vector_3  Vector_3;
  typedef typename boost::graph_traits<TriangulatedSurfaceMesh>::vertex_descriptor vertex_descriptor;

  typedef std::map<vertex_descriptor, FT> Vertex2FT_map;
  typedef boost::associative_property_map< Vertex2FT_map > Vertex2FT_property_map;

  typedef std::map<vertex_descriptor, Vector_3> Vertex2Vector_map;
  typedef boost::associative_property_map< Vertex2Vector_map > Vertex2Vector_property_map;
};


#endif
