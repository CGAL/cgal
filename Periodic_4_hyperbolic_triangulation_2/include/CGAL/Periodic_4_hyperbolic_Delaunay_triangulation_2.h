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
#include <map>

#include <CGAL/Square_root_2_field.h>
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_ds_vertex_base_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_ds_face_base_2.h>
#include <CGAL/Hyperbolic_word_4.h>


namespace CGAL {

template <  class GT,
            class TDS = Triangulation_data_structure_2<
              Triangulation_vertex_base_with_info_2<Hyperbolic_translation_info<std::string>, GT, Periodic_4_hyperbolic_triangulation_ds_vertex_base_2<GT> >,
              Triangulation_face_base_with_info_2<Periodic_4_hyperbolic_face_info_2, GT, Periodic_4_hyperbolic_triangulation_ds_face_base_2<> >
            >
          >
class Periodic_4_hyperbolic_Delaunay_triangulation_2 :
public Periodic_4_hyperbolic_triangulation_2<GT, TDS> {

  typedef Periodic_4_hyperbolic_Delaunay_triangulation_2<GT, TDS>   Self;
  typedef Periodic_4_hyperbolic_triangulation_2<GT, TDS>            Base;

public:
  typedef typename Base::Locate_type                         Locate_type;
  typedef typename Base::Geometric_traits                    Geometric_traits;
  typedef typename Base::Triangulation_data_structure        Triangulation_data_structure;
  typedef typename Base::Int                                 Int;
  typedef typename Base::Offset                              Offset;
  typedef typename Base::Circle_2                            Circle_2;
  typedef Circle_2                                           Circle;
  typedef typename Base::Point_2                             Point_2;
  typedef Point_2                                            Point;
  typedef typename Base::Segment_2                           Segment_2;
  typedef Segment_2                                          Segment;
  typedef typename Base::Triangle_2                          Triangle_2;
  typedef Triangle_2                                         Triangle;

  typedef typename Base::Periodic_point                      Periodic_point;
  typedef typename Base::Periodic_segment                    Periodic_segment;
  typedef typename Base::Periodic_triangle                   Periodic_triangle;  

  typedef typename Base::Vertex                              Vertex;
  typedef typename Base::Edge                                Edge;
  typedef typename Base::Face                                Face;

  typedef typename Base::Vertex_handle                       Vertex_handle;
  typedef typename Base::Face_handle                         Face_handle;

  typedef typename Base::size_type                           size_type;
  typedef typename Base::difference_type                     difference_type;

  typedef typename Base::Face_iterator                       Face_iterator;
  typedef typename Base::Edge_iterator                       Edge_iterator;
  typedef typename Base::Vertex_iterator                     Vertex_iterator;
  typedef typename Base::Face_circulator                     Face_circulator;
  typedef typename Base::Edge_circulator                     Edge_circulator;
  typedef typename Base::Vertex_circulator                   Vertex_circulator;

  typedef Face_iterator                                      All_faces_iterator;
  typedef Edge_iterator                                      All_edges_iterator;
  typedef Vertex_iterator                                    All_vertices_iterator;

  typedef Face_iterator                                      Finite_faces_iterator;
  typedef Edge_iterator                                      Finite_edges_iterator;
  typedef Vertex_iterator                                    Finite_vertices_iterator;

  typedef Periodic_4_hyperbolic_triangulation_unique_vertex_iterator_2<Self>
                                                             Unique_vertex_iterator;

private:
    typedef typename GT::FT                                  FT;
    typedef std::pair< Vertex_handle, Offset >               Virtual_vertex;
    typedef std::map<Vertex_handle, Virtual_vertex>          Virtual_vertex_map;
    typedef typename Virtual_vertex_map::const_iterator      Virtual_vertex_map_it;
    typedef std::map<Vertex_handle, std::vector<Vertex_handle > > 
                                                             Virtual_vertex_reverse_map;
    typedef typename Virtual_vertex_reverse_map::const_iterator
                                                             Virtual_vertex_reverse_map_it;
    typedef Triple< Vertex_handle, Vertex_handle, Vertex_handle >
                                                             Vertex_triple;

public:
    typedef Periodic_4_hyperbolic_triangulation_triangle_iterator_2<Self>
                                                             Periodic_triangle_iterator;
    typedef Periodic_4_hyperbolic_triangulation_segment_iterator_2<Self>
                                                             Periodic_segment_iterator;
    typedef Periodic_4_hyperbolic_triangulation_point_iterator_2<Self>
                                                             Periodic_point_iterator;

    typedef Point                                            value_type;
    typedef const value_type&                                const_reference;
    typedef Tag_false                                        Weighted_tag;

private:
  Geometric_traits              _gt;
  Triangulation_data_structure  _tds; 
  Circle                        _domain;


  void init_tds() {
    _tds.set_dimension(-2);
    v_offsets.reserve(14);
  }

protected:
  mutable std::vector<Vertex_handle> v_offsets;

public:

  Periodic_4_hyperbolic_Delaunay_triangulation_2(Geometric_traits gt) : 
    _gt(gt), _domain( Circle_2(Point(0,0), 1*1) ), _tds() {
    init_tds();
  }  

  Periodic_4_hyperbolic_Delaunay_triangulation_2(
    const Circle_2 domain = Circle_2(Point_2(FT(0),FT(0)), FT(1*1)), 
    const Geometric_traits &gt = Geometric_traits() ) :
    _gt(gt), _domain(domain), _tds() {
    //_gt.set_domain(_domain);
    init_tds();
  }

  Periodic_4_hyperbolic_Delaunay_triangulation_2(const Periodic_4_hyperbolic_Delaunay_triangulation_2& tr) :
    _gt(tr.geom_traits()),
    _domain(tr._domain) {
    CGAL_triangulation_expensive_postcondition(*this == tr);
  }


  Vertex_handle insert(const Point  &p,
                       Face_handle start = Face_handle() );

  Vertex_handle insert(const Point& p,
                       Locate_type lt,
                       Face_handle loc, int li );

  template < class InputIterator >
  std::ptrdiff_t insert(InputIterator first, 
                        InputIterator last, 
                        bool is_large_point_set = true);


};  // class Periodic_4_hyperbolic_Delaunay_triangulation_2

    
template < class Gt, class Tds >
inline
typename Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
insert(const Point  &p,  Face_handle start)
{
  std::cout << "Inseritng now! Yeah!" << std::endl;
  return Vertex_handle();
}


template < class Gt, class Tds >
inline
typename Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
insert(const Point  &p, Locate_type lt, Face_handle loc, int li)
{
  std::cout << "Inseritng now! Right!" << std::endl;
  return Vertex_handle();
}


} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
