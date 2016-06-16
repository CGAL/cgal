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
              Periodic_4_hyperbolic_triangulation_ds_vertex_base_2<GT>,
              Periodic_4_hyperbolic_triangulation_ds_face_base_2<GT> 
            >
          >
class Periodic_4_hyperbolic_Delaunay_triangulation_2 :
public Periodic_4_hyperbolic_triangulation_2<GT, TDS> {

  typedef Periodic_4_hyperbolic_Delaunay_triangulation_2<GT, TDS>   Self;
  typedef Periodic_4_hyperbolic_triangulation_2<GT, TDS>            Base;

public:

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2  
  using Base::cw;
  using Base::ccw;
  using Base::geom_traits;
#endif

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
  typedef typename Base::Line_face_circulator                Line_face_circulator;

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

  bool locally_Delaunay(const Face_handle&, int);
  void propagating_flip(Face_handle&, int);
  void restore_Delaunay(Vertex_handle);
  void flip_single_edge(Face_handle, int);
  bool flippable(Face_handle, int);

};  // class Periodic_4_hyperbolic_Delaunay_triangulation_2




template<class Gt, class Tds>
bool Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::flippable(Face_handle f, int i)
{
  Face_handle nb = f->neighbor(i);
  int j = nb->index(f);

  const Point *p[4];

  p[0] = &f->vertex(i)->point();      // i
  p[1] = &nb->vertex(j)->point();     // opposite
  p[2] = &f->vertex(ccw(i))->point(); // ccw
  p[3] = &f->vertex(cw(i))->point();  // cw

  if (f->has_zero_offsets() && nb->has_zero_offsets()) {
      if (orientation(*p[0], *p[1], *p[2]) == LEFT_TURN)
        return false;
      if (orientation(*p[0], *p[1], *p[3]) == RIGHT_TURN)
        return false;
  } else {
      Offset off[4];
      off[0] = f->offset(i);
      off[1] = f->neighbor_offset(j).append(nb->offset(j)); 
      off[2] = f->offset(ccw(i));
      off[3] = f->offset(cw(i));

      if (orientation(*p[0], *p[1], *p[2], off[0], off[1], off[2]) == LEFT_TURN)
        return false;
      if (orientation(*p[0], *p[1], *p[3], off[0], off[1], off[3]) == RIGHT_TURN)
        return false;
  }

  return true;
}



template<class Gt, class Tds>
void Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::flip_single_edge(Face_handle f, int i)
{
  CGAL_triangulation_precondition(f != Face_handle());
  CGAL_triangulation_precondition(i == 0 || i == 1 || i == 2);
  CGAL_triangulation_precondition(dimension() == 2);

  CGAL_triangulation_precondition(flippable(f, i));

  if (f->has_zero_offsets() && f->neighbor(i)->has_zero_offsets()) {
    _tds.flip(f, i);
    return;
  } else {
      Face_handle nb = f->neighbor(i);
      int nb_idx = nb->index(f);
      
      Vertex_handle vh[] = { f->vertex(i),
                             f->vertex(ccw(i)),
                            nb->vertex(nb_idx),
                             f->vertex(cw(i))    };
      Offset o[]         = { f->offset(i),
                             f->offset(ccw(i)),
                            nb->offset(nb_idx),
                             f->offset(cw(i))    };
      Offset no[]        = { f->neighbor_offset(cw(i)),
                            nb->neighbor_offset(ccw(nb_idx)),
                            nb->neighbor_offset(cw(nb_idx)),
                             f->neighbor_offset(ccw(i))  };

      _tds.flip(f, i);
      
      nb = f->neighbor(ccw(i));
      nb_idx = nb->index(f);

      f->set_offsets();
      nb->set_offsets();

      f->set_offset(i, o[0]);
      f->set_neighbor_face_offset(i, no[1]);
      f->set_offset(ccw(i), o[1]);
      f->set_neighbor_face_offset(ccw(i), Offset());
      f->set_offset(cw(i), o[2]); 
      f->set_neighbor_face_offset(cw(i), no[0]);

      nb->set_offset(nb_idx, o[3]);
      nb->set_neighbor_face_offset(nb_idx, Offset());
      nb->set_offset(ccw(nb_idx), o[0]);
      nb->set_neighbor_face_offset(ccw(nb_idx), no[2]);
      nb->set_offset(cw(nb_idx), o[2]);
      nb->set_neighbor_face_offset(cw(nb_idx), no[3]);
  }


}


    


template < class Gt, class Tds >
void
Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
restore_Delaunay(Vertex_handle v)
{
  Face_handle f = v->face();
  Face_handle next;
  int i;
  Face_handle start(f);
  do
    {
      i = f->index(v);
      next = f->neighbor(ccw(i));  // turn ccw around v
      propagating_flip(f, i);
      f = next;
    }
  while(next != start);
}


template < class Gt, class Tds >
void
Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
propagating_flip(Face_handle& f, int i)
{
  Face_handle nb = f->neighbor(i);

  if (locally_Delaunay(f, i))
    return;

  this->flip_single_edge(f, i);
  propagating_flip(f, i);
  i = nb->index(f->vertex(i));
  propagating_flip(nb, i);
}



template < class Gt, class Tds >
bool
Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
locally_Delaunay(const Face_handle &f, int nbi)
{
  CGAL_BRANCH_PROFILER("locally_Delaunay(), simplicity check failures", tmp);

  Face_handle nb = f->neighbor(nbi);

  bool simplicity_criterion = f->has_zero_offsets() && nb->has_zero_offsets();

  const Point *p[4];
  for (int index = 0; index < 3; ++index)
    {
      p[index]   = &nb->vertex(index)->point();
    }
  p[3]   = &f->vertex(nbi)->point();

  Oriented_side os;
  if (simplicity_criterion)
    {
      // No periodic offsets
      os = side_of_oriented_circle(*p[0], *p[1], *p[2], *p[3]);
    }
  else
    {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);

      Offset off[4];

      for (int index = 0; index < 3; ++index)
        {
          off[index] = nb->offset(index);
        }
      off[3] = f->neighbor_offset(nbi);

      os = side_of_oriented_circle( off[0].apply(*p[0]), 
                                    off[1].apply(*p[1]), 
                                    off[2].apply(*p[2]), 
                                    off[3].apply(*p[3]) );
    }

  return (ON_POSITIVE_SIDE != os);
}






template < class Gt, class Tds >
inline
typename Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
insert(const Point  &p,  Face_handle start)
{


  typedef typename Gt::Side_of_fundamental_octagon Side_of_fundamental_octagon;
  
  Side_of_fundamental_octagon check = Side_of_fundamental_octagon();
  CGAL::Bounded_side side = check(p);

  if (side != CGAL::ON_UNBOUNDED_SIDE) {

    if ( start == Face_handle() ) {
      Locate_type lt;
      int li;
      start = this->euclidean_visibility_locate(p, lt, li);
    }

    cout << "Point inserted in face " << start->get_number() << ", vertices: " << endl;
    for (int i = 0; i < 3; i++) 
      cout << start->vertex(i)->idx() << " with offset " << start->offset(i) << ", " << endl;
    cout << endl;
    cout << "Neighbor offsets: " << start->neighbor_offset(0) << ", " << start->neighbor_offset(1) << ", " << start->neighbor_offset(2) << endl;

    int next_number = this->number_of_faces();
    Vertex_handle ov[] = {  start->vertex(0),
                            start->vertex(1),
                            start->vertex(2) };
    Offset        oo[] = {  start->offset(0),
                            start->offset(1),
                            start->offset(2) };
    Offset       ono[] = {  start->neighbor_offset(0),
                            start->neighbor_offset(1),
                            start->neighbor_offset(2) };

    Vertex_handle new_vertex = this->insert_in_face(p, start);
    new_vertex->set_point( p );
    new_vertex->set_idx( this->number_of_vertices() );
    cout << "New vertex has id " << new_vertex->idx() << endl;

    // Assign offsets
    Face_circulator iface = _tds.incident_faces(new_vertex), end(iface);
    if (iface != 0) {
      do {
        
        // Null all offsets
        iface->set_offsets();

        if (iface->get_number() == -1) {
          iface->set_number(next_number++);
        }
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            if (iface->vertex(i) == ov[j]) {
              iface->set_offset(i, oo[j]);
              if (iface->vertex(Triangulation_cw_ccw_2::cw(i)) == ov[Triangulation_cw_ccw_2::cw(j)]) {
                iface->set_neighbor_face_offset(Triangulation_cw_ccw_2::ccw(i), ono[Triangulation_cw_ccw_2::ccw(j)]);
              } 
              else if (iface->vertex(Triangulation_cw_ccw_2::ccw(i)) == ov[Triangulation_cw_ccw_2::ccw(j)]) {
                iface->set_neighbor_face_offset(Triangulation_cw_ccw_2::cw(i), ono[Triangulation_cw_ccw_2::cw(j)]);
              }
            }
          }
        }

      } while (++iface != end);
    
      cout << "New faces created: " << endl << "------------------------" << endl;
      do {
        cout << "Face " << iface->get_number() << ", vertices: ";
        for (int i = 0; i < 3; i++) {
          cout << iface->vertex(i)->idx() << " with offset " << iface->offset(i) << ", ";
        }
        cout << endl;
        cout << "Neighbor offsets: " << iface->neighbor_offset(0) << ", " << iface->neighbor_offset(1) << ", " << iface->neighbor_offset(2) << endl;
      } while (++iface != end);
      cout << "--------- end ----------" << endl << endl;
    }

    restore_Delaunay(new_vertex);

    return new_vertex;
  } else {
    return Vertex_handle();
  }

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
