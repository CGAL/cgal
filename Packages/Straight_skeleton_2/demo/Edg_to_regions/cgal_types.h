// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Radu Ursu

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <boost/shared_ptr.hpp>

namespace demo
{

typedef double                     NT;
typedef CGAL::Cartesian<NT>        R;
typedef CGAL::Filtered_kernel<R>   K;
typedef K::Point_2                 Point;
typedef CGAL::Segment_2<K>         Segment;
typedef CGAL::Polygon_2<K>         Polygon;
typedef boost::shared_ptr<Polygon> PolygonPtr;
typedef std::vector<PolygonPtr>    Region ;
typedef boost::shared_ptr<Region>  RegionPtr;
typedef std::vector<RegionPtr>     RegionList ;

typedef CGAL::Triangulation_data_structure_2<CGAL::Triangulation_vertex_base_2<K>
                                            ,CGAL::Constrained_triangulation_face_base_2<K> 
                                            > TDS ;
                                            
typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS,CGAL::Exact_predicates_tag> CDT;

typedef CDT::Vertex_handle Vertex_handle ;
typedef CDT::Face_handle   Face_handle ;

}
