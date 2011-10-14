// Copyright (c) 2010 CNRS, LIRIS, http://liris.cnrs.fr/, All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Cell_attribute.h>
#include <iostream>
#include <algorithm>

struct Average_functor
{
  template<class CellAttribute>
  void operator()(CellAttribute& ca1,const CellAttribute& ca2)
  { ca1.info()=(ca1.info()+ ca2.info())/2; }
};

struct Myitem
{
  template<class Refs>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<3, Refs > Dart;
    
    typedef CGAL::Cell_attribute_with_point< Refs, int, CGAL::Tag_true, 
					     Average_functor >
    Vertex_attribute;
    
    typedef CGAL::cpp0x::tuple<Vertex_attribute> Attributes;
  };
};

typedef CGAL::Linear_cell_complex_traits<3,
           CGAL::Exact_predicates_inexact_constructions_kernel> Traits;
typedef CGAL::Combinatorial_map_with_points<3,3,Traits,Myitem> LCC_3;
typedef LCC_3::Dart_handle                             Dart_handle;
typedef LCC_3::Point                                   Point;
typedef LCC_3::FT                                      FT;

Dart_handle make_iso_cuboid(LCC_3& lcc, const Point& basepoint, FT lg)
{
	return make_hexahedron(lcc,basepoint,
												 LCC_3::Construct_translated_point()(basepoint,LCC_3::Vector(lg,0,0)),
												 LCC_3::Construct_translated_point()(basepoint,LCC_3::Vector(lg,lg,0)),
												 LCC_3::Construct_translated_point()(basepoint,LCC_3::Vector(0,lg,0)),
												 LCC_3::Construct_translated_point()(basepoint,LCC_3::Vector(0,lg,lg)),
												 LCC_3::Construct_translated_point()(basepoint,LCC_3::Vector(0,0,lg)),
												 LCC_3::Construct_translated_point()(basepoint,LCC_3::Vector(lg,0,lg)),
												 LCC_3::Construct_translated_point()(basepoint,LCC_3::Vector(lg,lg,lg)));
}

int main()
{
  LCC_3 lcc;
  
  // Create 2 cubes.
  Dart_handle d1 = make_iso_cuboid(lcc, Point(-2, 0, 0), 1);
  Dart_handle d2 = make_iso_cuboid(lcc, Point(0, 0, 0), 1);

  // Set the color of all vertices of the first cube to 1
  for (LCC_3::One_dart_per_incident_cell_range<0, 3>::iterator 
	 it=lcc.one_dart_per_incident_cell<0,3>(d1).begin(), 
	 itend=lcc.one_dart_per_incident_cell<0,3>(d1).end(); it!=itend; ++it)
    { LCC_3::vertex_attribute(it)->info()=1; }
  
  // Set the color of all vertices of the second cube to 19
  for (LCC_3::One_dart_per_incident_cell_range<0, 3>::iterator it=
	 lcc.one_dart_per_incident_cell<0,3>(d2).begin(),
	 itend=lcc.one_dart_per_incident_cell<0,3>(d2).end(); it!=itend; ++it)
    { LCC_3::vertex_attribute(it)->info()=19; }
  
  // 3-Sew the two cubes along one facet
  lcc.sew<3>(d1->beta(1)->beta(1)->beta(2), d2->beta(2));

  // Barycentric triangulation of the facet between the two cubes.
  Dart_handle d3=insert_center_cell_0_in_cell_2(lcc, d2->beta(2));

  // Set the color of the new vertex to 5.
  LCC_3::vertex_attribute(d3)->info()=5;
  
  // Display all the vertices of the map.
  for (LCC_3::One_dart_per_cell_range<0>::iterator 
	 it=lcc.one_dart_per_cell<0>().begin(), 
	 itend=lcc.one_dart_per_cell<0>().end(); 
       it!=itend; ++it)
    {
      std::cout<<"point: "<<LCC_3::point(it)<<", "
               <<"color: "<<LCC_3::vertex_attribute(it)->info()
               <<std::endl;
    }

  return EXIT_SUCCESS;
}
