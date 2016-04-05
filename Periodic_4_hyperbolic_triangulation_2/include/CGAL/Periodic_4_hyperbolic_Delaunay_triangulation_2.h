// Copyright (c) 2016  INRIA Nancy - Grand Est (France).
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
// $URL: $
// $Id:  $
//
//
// Author(s)     : Iordan Iordanov

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H

#include <vector>

#include <CGAL/Square_root_2_field.h>
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>
#include <CGAL/Delaunay_hyperbolic_triangulation_2.h>

namespace CGAL {
    
template <  class Gt, 
            class Tds = Triangulation_data_structure_2 <
                        Triangulation_vertex_base_2<Gt>, 
                        Triangulation_face_base_with_info_2<Hyperbolic_face_info_2, Gt> > >
class Periodic_4_hyperbolic_Delaunay_triangulation_2 : public Delaunay_hyperbolic_triangulation_2<Gt, Tds> {

public:
    typedef Delaunay_hyperbolic_triangulation_2<Gt, Tds>                Base;
    typedef Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>     Self;

    typedef typename Base::Face_base                                    Face_base;
    typedef typename Base::size_type                                    size_type;
    typedef typename Base::Vertex_handle                                Vertex_handle;
    typedef typename Base::Face_handle                                  Face_handle;
    typedef typename Base::Edge                                         Edge;

    typedef typename Base::Edge_circulator                              Edge_circulator;
    typedef typename Base::Face_circulator                              Face_circulator;
    typedef typename Base::Vertex_circulator                            Vertex_circulator;
  
    typedef typename Base::All_vertices_iterator                        All_vertices_iterator;
    typedef typename Base::All_edges_iterator                           All_edges_iterator;
    typedef typename Base::All_faces_iterator                           All_faces_iterator;
 
    typedef Gt                                                          Geom_traits;
    typedef typename Geom_traits::FT                                    FT;
    typedef typename Geom_traits::Point_2                               Point;
    typedef typename Geom_traits::Segment_2                             Segment;

    // Periodic stuff
    typedef Hyperbolic_octagon_translation_matrix                       Offset;
    typedef std::pair<Point, Offset>                                    Periodic_point;
    typedef array< std::pair<Point,Offset>, 2>                          Periodic_segment;
    typedef array< std::pair<Point,Offset>, 3>                          Periodic_triangle;


    Periodic_4_hyperbolic_Delaunay_triangulation_2(const Gt& gt = Gt())
    : Delaunay_hyperbolic_triangulation_2<Gt,Tds>(gt) {}
  
    Periodic_4_hyperbolic_Delaunay_triangulation_2(
         const Delaunay_hyperbolic_triangulation_2<Gt,Tds> &tr)
    : Delaunay_hyperbolic_triangulation_2<Gt,Tds>(tr) {   
        CGAL_triangulation_postcondition( this->is_valid() );
    }

    void insert_dummy_points(std::vector<Point>&);

};

    
} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
