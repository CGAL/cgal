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
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Circular_kernel_intersections.h>

#include <iterator>
#include <CGAL/intersections.h>
#include <CGAL/result_of.h>
#include <CGAL/iterator.h>
#include <boost/bind.hpp>



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


    class Dummy_point {
    private:
      Point _pt;
      bool  _is_inserted;
      Vertex_handle _vh;

    public:

      Dummy_point(FT x, FT y): _pt(x,y), _is_inserted(true) {}
      Dummy_point(Point p):    _pt(p),   _is_inserted(true) {}

      Point operator()()      {  return _pt;          }
      bool  is_inserted()     {  return _is_inserted; }
      Vertex_handle vertex()  {  return _vh;          }
      void  set_inserted(bool val)      { _is_inserted = val; }
      void  set(Point val)              { _pt = val;          }
      void  set_vertex(Vertex_handle v) { _vh = v;            }
    };

    std::vector<Dummy_point> dummy_points;

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
    int f_cnt, v_cnt, n_dpt;


  private:
    bool is_removable(Vertex_handle v, Delaunay_triangulation_2<GT, TDS>& dt, std::map<Vertex_handle, Vertex_handle>& vmap);

  public:
    Periodic_4_hyperbolic_Delaunay_triangulation_2(Geometric_traits gt) : 
    Periodic_4_hyperbolic_triangulation_2<GT, TDS>(gt) { n_dpt = 14; }  

    Periodic_4_hyperbolic_Delaunay_triangulation_2(
      const Circle_2 domain = Circle_2(Point_2(FT(0),FT(0)), FT(1*1)), 
      const Geometric_traits &gt = Geometric_traits() ) :
    Periodic_4_hyperbolic_triangulation_2<GT, TDS>(domain, gt) { n_dpt = 14; }

    Periodic_4_hyperbolic_Delaunay_triangulation_2(const Periodic_4_hyperbolic_Delaunay_triangulation_2& tr) :
    Periodic_4_hyperbolic_triangulation_2<GT, TDS>(tr) { }

    std::vector<Vertex_handle> insert_dummy_points(bool rational = true);

    Vertex_handle insert(const Point  &p,
     Face_handle start = Face_handle() );

    template < class InputIterator >
    std::ptrdiff_t insert(InputIterator first, 
      InputIterator last);

    void remove(Vertex_handle v);

    int number_of_dummy_points() { return n_dpt; }

    bool _side_of_octagon( const Face_handle& fh, const Offset& offset) const {
      int cnt = 0;
      typename GT::Side_of_fundamental_octagon side;
      for (int j = 0; j < 3; j++) {
        Offset o = offset.inverse().append(fh->vertex(j)->get_offset());
        Point  p = o.apply( fh->vertex(j)->point() );
        if ( side(p) == CGAL::ON_UNBOUNDED_SIDE ) {
          if ( p.y() + tan(CGAL_PI / FT(8))*p.x() > 0 ) {
            cnt++;
          } else {
          }
        }
      }
      return (cnt == 0);
    }

};  // class Periodic_4_hyperbolic_Delaunay_triangulation_2


template <class Gt>
typename Gt::FT dist(typename Gt::Point_2 p, typename Gt::Point_2 q) {
  typename Gt::FT r = (p.x() - q.x())*(p.x() - q.x()) + (p.y() - q.y())*(p.y() - q.y());
  return sqrt(r);
}


template <class Gt>
typename Gt::FT 
hyperbolic_diameter(typename Gt::Circle_2 c) {
  typedef typename Gt::FT           FT;
  typedef typename Gt::Kernel       K;
  typedef typename Gt::Point_2      Point;
  typedef typename Gt::Line_2       Line;
  typedef typename Gt::Circle_2     Circle;


  typedef CGAL::Exact_circular_kernel_2        CircK;
  typedef CGAL::Point_2<CircK>                 Pt2;
  typedef CGAL::Circle_2<CircK>                Circ2;
  typedef CGAL::Line_2<CircK>                  Line2;
  typedef std::pair<CGAL::Circular_arc_point_2<CircK>, unsigned> IsectOutput;

  typedef CGAL::Dispatch_output_iterator< 
                CGAL::cpp11::tuple<IsectOutput>, 
                CGAL::cpp0x::tuple< std::back_insert_iterator<std::vector<IsectOutput> > > > Dispatcher;

  std::vector<IsectOutput> res0, res1;
  Dispatcher disp1 = CGAL::dispatch_output<IsectOutput>( std::back_inserter(res1) ); 
  Dispatcher disp0 = CGAL::dispatch_output<IsectOutput>( std::back_inserter(res0) ); 

  Pt2 ct(to_double(c.center().x()), to_double(c.center().y()));
  Pt2 p0(0, 0);
  double r = to_double(c.squared_radius());

  Line2 ell(p0, ct);
  Circ2 c2(ct, r);
  Circ2 c0(p0, 1);

  if (ell.is_degenerate()) {
    //cout << "\tThis is degenerate case!" << endl;
    return 5.;
  } else {
    CGAL::intersection(c0, ell, disp0);
    CGAL::intersection(c2, ell, disp1);
  }
  Point a(to_double(res0[0].first.x()), to_double(res0[0].first.y()));
  Point b(to_double(res0[1].first.x()), to_double(res0[1].first.y()));

  Point p(to_double(res1[0].first.x()), to_double(res1[0].first.y()));
  Point q(to_double(res1[1].first.x()), to_double(res1[1].first.y()));
  
  FT aq = dist<Gt>(a, q);
  FT pb = dist<Gt>(p, b);
  FT ap = dist<Gt>(a, p);
  FT qb = dist<Gt>(q, b);

  //cout << "aq = " << aq << ", pb = " << pb << " | ap = " << ap << ", qb = " << qb << endl;

  double hyperdist = fabs(log(to_double((aq*pb)/(ap*qb))));

  return hyperdist;
}


template <class Gt, class Tds>
bool
Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt,Tds>::
is_removable(Vertex_handle v, Delaunay_triangulation_2<Gt,Tds>& dt, std::map<Vertex_handle, Vertex_handle>& vmap) {
  
  typedef typename Gt::FT                             FT;
  typedef typename Gt::Circle_2                       Circle;
  typedef Delaunay_triangulation_2<Gt, Tds>           Delaunay;
  typedef typename Delaunay::Finite_faces_iterator    Finite_Delaunay_faces_iterator;

  // This is a rational approximation of the limit (1/4 of the squared systole)
  Circle lc(dummy_points[0](),
            dummy_points[13](),
            Offset(1,6,3).apply(dummy_points[13]()));
  FT lim = hyperbolic_diameter<Gt>(lc)/FT(2);  

  // This is the exact value of the limit (1/4 of the squared systole)
  //FT lim = FT( pow(acosh(1. + sqrt(2.)),2.)/4. );

  std::vector<Vertex_handle> bdry_verts;
  Face_circulator nbf(tds().incident_faces(v)), done(nbf);
  do {
    int idx = nbf->index(v);
    Offset off = nbf->offset(idx).inverse();
    off = off*nbf->offset(ccw(idx));
    Vertex_handle thisv = nbf->vertex(ccw(idx));
    bdry_verts.push_back(thisv);
    Point pt = off.apply(thisv->point());
    Vertex_handle new_v = dt.insert(pt);
    vmap.insert(std::pair<Vertex_handle, Vertex_handle>(new_v, thisv));
  } while (++nbf != done);


  int n_verts = bdry_verts.size();
  FT max_diam = FT(0);
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
      Circle c(fit->vertex(0)->point(), 
               fit->vertex(1)->point(), 
               fit->vertex(2)->point());
      FT diam = hyperbolic_diameter<Gt>(c);
      if (max_diam < diam) {
        max_diam = diam;
      }
    }
  }

  if (max_diam < lim) {
    return true;
  } else {
    return false;
  }
}


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

    for (int i = 0; i < dummy_points.size(); i++) {
      if (dummy_points[i].is_inserted()) {
        typedef Delaunay_triangulation_2<Gt, Tds>           Delaunay;
        Delaunay dt;
        std::map<Vertex_handle, Vertex_handle> vmap;

        if (is_removable(dummy_points[i].vertex(), dt, vmap)) {
          //cout << "Removing dummy point " << i << endl;
          remove(dummy_points[i].vertex());
          dummy_points[i].set_inserted(false);
        }
      }
    }

    n_dpt = 0;
    for (int i = 0; i < dummy_points.size(); i++) {
      if (dummy_points[i].is_inserted())
        n_dpt++;
    }

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
  Delaunay dt;
  std::map<Vertex_handle, Vertex_handle> vmap;

  if (is_removable(v, dt, vmap)) {
  
    typedef typename Delaunay::Finite_faces_iterator    Finite_Delaunay_faces_iterator;
    typedef std::pair<Face_handle, int>                 Neighbor_pair;
    typedef std::pair<Edge, Neighbor_pair>              Edge_neighbor;


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
      bdry_verts.push_back(nb->vertex(ccw(idx)));

      nb->store_offsets(nb->offset(idx).inverse());
      nbrs.push_back(nb);
      nb++;
    } while(nb != done);

    for (int i = 0; i < bdry_edges.size(); i++) {
      Edge e = bdry_edges[i];
      Face_handle f = e.first;
      int j = e.second;
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

  } else { // is not removable
    //cout << "Vertex " << v->idx() << " cannot be removed!" << endl;
  }

}


} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
