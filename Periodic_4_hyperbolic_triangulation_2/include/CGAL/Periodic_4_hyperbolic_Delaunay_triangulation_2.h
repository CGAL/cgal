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
#include <CGAL/Hyperbolic_octagon_word_4.h>

#include <CGAL/Delaunay_triangulation_2.h>

namespace CGAL {

  template <  class GT,
  class TDS = Triangulation_data_structure_2<
  Periodic_4_hyperbolic_triangulation_ds_vertex_base_2<GT>,
  Periodic_4_hyperbolic_triangulation_ds_face_base_2<GT> 
  >
  >
  class Periodic_4_hyperbolic_Delaunay_triangulation_2: public Periodic_4_hyperbolic_triangulation_2<GT, TDS> {

    typedef Periodic_4_hyperbolic_Delaunay_triangulation_2<GT, TDS>   Self;
    typedef Periodic_4_hyperbolic_triangulation_2<GT, TDS>            Base;

  public:

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2  
    using Base::cw;
    using Base::ccw;
    using Base::geom_traits;
    using Base::tds;
#endif

    typedef typename Base::Locate_type                         Locate_type;
    typedef typename Base::Geometric_traits                    Geometric_traits;
    typedef typename Base::Triangulation_data_structure        Triangulation_data_structure;
    //typedef typename Base::Int                                 Int;
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

  protected:
    mutable std::vector<Vertex_handle> v_offsets;

        int f_cnt, v_cnt;

  public:

    Periodic_4_hyperbolic_Delaunay_triangulation_2(Geometric_traits gt) : 
    Periodic_4_hyperbolic_triangulation_2<GT, TDS>(gt) { }  

    Periodic_4_hyperbolic_Delaunay_triangulation_2(
      const Circle_2 domain = Circle_2(Point_2(FT(0),FT(0)), FT(1*1)), 
      const Geometric_traits &gt = Geometric_traits() ) :
    Periodic_4_hyperbolic_triangulation_2<GT, TDS>(domain, gt) { }

    Periodic_4_hyperbolic_Delaunay_triangulation_2(const Periodic_4_hyperbolic_Delaunay_triangulation_2& tr) :
    Periodic_4_hyperbolic_triangulation_2<GT, TDS>(tr) { }

    std::vector<Vertex_handle> insert_dummy_points(bool rational = true);

    Vertex_handle insert(const Point  &p,
     Face_handle start = Face_handle() );

    template < class InputIterator >
    std::ptrdiff_t insert(InputIterator first, 
      InputIterator last);

    void remove(Vertex_handle v);

    bool _side_of_octagon( const Face_handle& fh, const Offset& offset) const {
      int cnt = 0;
      typename GT::Side_of_fundamental_octagon side;
      for (int j = 0; j < 3; j++) {
        Offset o = offset.inverse().append(fh->vertex(j)->get_offset());
        Point  p = o.apply( fh->vertex(j)->point() );
        if ( side(p) == CGAL::ON_UNBOUNDED_SIDE ) {
          if ( p.y() + (CGAL_PI / FT(8))*p.x() > 0 ) {
            cnt++;
          } else {
          }
        }
      }
      return (cnt == 0);
    }

};  // class Periodic_4_hyperbolic_Delaunay_triangulation_2




template < class Gt, class Tds >
inline
typename Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
insert(const Point  &p,  Face_handle start) {

  Vertex_handle v;

  typedef typename Gt::Side_of_fundamental_octagon Side_of_fundamental_octagon;

  Side_of_fundamental_octagon check = Side_of_fundamental_octagon();
  CGAL::Bounded_side side = check(p);

  if (side != CGAL::ON_UNBOUNDED_SIDE) {
      Offset loff;
      if ( start == Face_handle() ) {
      Locate_type lt;
      int li;
      start = this->euclidean_visibility_locate(p, lt, li, loff);
      if (lt == Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::VERTEX) {
        return Vertex_handle();
      }
    }

    std::vector<Face_handle> faces;
    this->find_conflicts(start, p, loff, std::back_inserter(faces));
    v = this->insert_in_hole(p, faces.begin(), faces.end());
    v->set_idx(v_cnt++);

    Face_circulator ifc = tds().incident_faces(v), done(ifc);
    do {
      ifc->restore_offsets(loff);
      ifc->tds_data().clear();
      ifc->make_canonical();
      ifc->set_number(f_cnt++);
    } while (++ifc != done);

    Vertex_circulator ivc = tds().incident_vertices(v), done_v(ivc);
    do {
      ivc->remove_offset();
    } while (++ivc != done_v);

    CGAL_triangulation_assertion(this->is_valid());

    return v;
  }

  return Vertex_handle();
}


//------------------------------------------------------

template < class Gt, class Tds >
void
Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
remove(Vertex_handle v) {

  typedef Delaunay_triangulation_2<Gt, Tds>           Delaunay;
  typedef typename Delaunay::Finite_faces_iterator    Finite_Delaunay_faces_iterator;
  typedef std::pair<Face_handle, int>                 Neighbor_pair;
  typedef std::pair<Edge, Neighbor_pair>              Edge_neighbor;


  Delaunay dt;

  std::vector<Edge>             bdry_edges;
  std::vector<Vertex_handle>    bdry_verts;
  std::map<Edge, Neighbor_pair> bdry_nbrs;

  Face_circulator nb = tds().incident_faces(v), done(nb);
  std::vector<Face_handle> nbrs;
  do {
    int idx = nb->index(v);
    Edge e = Edge(nb, idx);
    bdry_edges.push_back(e);
    Face_handle nbf = nb->neighbor(idx);
    int nidx = 0;
    if (nbf->neighbor(1) == nb) nidx = 1;
    if (nbf->neighbor(2) == nb) nidx = 2;
    CGAL_triangulation_assertion(nbf->neighbor(nidx) == nb);
    bdry_nbrs.insert(Edge_neighbor(e, Neighbor_pair(nbf, nidx)));

    nb->store_offsets(nb->offset(idx).inverse());
    nbrs.push_back(nb);
    nb++;
  } while(nb != done);

  for (int i = 0; i < bdry_edges.size(); i++) {
    Edge e = bdry_edges[i];
    Face_handle f = e.first;
    int j = e.second;
  }

  std::map<Vertex_handle, Vertex_handle> vmap;

  for (int i = 0; i < bdry_edges.size(); i++) {
    Face_handle f = bdry_edges[i].first;
    int j = bdry_edges[i].second;
    Vertex_handle vh = f->vertex(ccw(j));
    bdry_verts.push_back(vh);
    Vertex_handle nv = dt.insert(vh->get_offset().apply(vh->point()));
    vmap.insert(std::pair<Vertex_handle, Vertex_handle>(nv, vh));
  }

  int n_verts = bdry_verts.size();
  std::vector<Face_handle> new_f;
  for (Finite_Delaunay_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); fit++) {
    bool is_good = true;
    for (int i = 0; i < 3; i++) {
      Vertex_handle this_v = vmap[fit->vertex(i)];
      Vertex_handle prev_v = bdry_verts[n_verts - 1];
      Vertex_handle curr_v = bdry_verts[0];
      for (int j = 1; curr_v != this_v; j = (j+1)%n_verts) {
        prev_v = curr_v;
        curr_v = bdry_verts[j];
      }
      if (vmap[fit->vertex(ccw(i))] == prev_v) {
        is_good = false;
        break;
      }
    }

    if (is_good) {
      Face_handle f = tds().create_face();
      f->set_number(f_cnt++);
      for (int j = 0; j < 3; j++) {
        f->set_vertex(j, vmap[fit->vertex(j)]);
      }
      new_f.push_back(f);
    }
  }


  int internb = 0;
  int bdrynb = 0;
  for (int i = 0; i < new_f.size(); i++) {
    for (int k = 0; k < 3; k++) {
      bool found_bdry = false;
      for (int j = 0; j < bdry_verts.size(); j++) { 
        if (new_f[i]->vertex(ccw(k)) == bdry_verts[j] && 
            new_f[i]->vertex(cw(k)) == bdry_verts[(j+1)%n_verts]) {
          found_bdry = true;
          Neighbor_pair nb = bdry_nbrs[bdry_edges[j]];
          Face_handle nbf = nb.first;
          int nbidx = nb.second;
          tds().set_adjacency(nbf, nbidx, new_f[i], k);
          bdrynb++;
          break;
        } 
      }
      if (!found_bdry) {
        for (int l = 0; l < new_f.size(); l++) {
          if (l == i) continue;
          for (int j = 0; j < 3; j++) {
            if (new_f[i]->vertex(ccw(k)) == new_f[l]->vertex(cw(j)) &&
                new_f[i]->vertex(cw(k))  == new_f[l]->vertex(ccw(j)) ) {
              tds().set_adjacency(new_f[i], k, new_f[l], j);
              internb++;
              break;
            }
          }
        }
      }
    }
  }


  for (int j = 0; j < new_f.size(); j++) {
    for (int i = 0; i < 3; i++) {
      new_f[j]->vertex(i)->set_face(new_f[j]);
    }
    new_f[j]->restore_offsets();
    new_f[j]->make_canonical();
  }

  for (int j = 0; j < bdry_edges.size(); j++) {
    Face_handle f = bdry_edges[j].first;
    int i = bdry_edges[j].second;
    f->vertex(ccw(i))->remove_offset();
  }

  for (int i = 0; i < nbrs.size(); i++) {
    tds().delete_face(nbrs[i]);
  }
  tds().delete_vertex(v);

  CGAL_triangulation_assertion(this->is_valid());
}


} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
