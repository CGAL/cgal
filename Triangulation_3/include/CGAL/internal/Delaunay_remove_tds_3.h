// Copyright (c) 2001  INRIA Sophia-Antipolis (France).
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
// Author(s)     :  Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//                  Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

#ifndef CGAL_INTERNAL_DELAUNAY_REMOVE_TDS_3_H
#define CGAL_INTERNAL_DELAUNAY_REMOVE_TDS_3_H

#include <CGAL/basic.h>

#include <map>
#include <CGAL/Triangulation_ds_face_base_2.h>
#include <CGAL/Triangulation_ds_vertex_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>

namespace CGAL { namespace internal {

// Should this be documented as a Triangulation_ds_vertex_base_with_info_2 ?
template < class I, class Vb = Triangulation_ds_vertex_base_2<> >
class Delaunay_remove_tds_vertex_3_2
  : public Vb
{
  I _info;
public :

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other   Vb2;
    typedef Delaunay_remove_tds_vertex_3_2<I, Vb2>          Other;
  };

  void set_info(const I &info) { _info = info; }

  const I & info() const { return _info; }
};


/* We derive the face class, because we need additional pointers to
   a successor and predecessor.

   We want to change the order to avoid looking at faces too often.
   When we make an operation (flip of an edge / removal of a vertex)
   we mark the adjacent four / three edges and move them right
   after the current face of the traversal. As they are marked
   they must be considered later, and no faces with unmarked edges
   are traversed to reach them.
*/

template < class I, class Fb = Triangulation_ds_face_base_2<> >
class Delaunay_remove_tds_face_3_2
  : public Fb
{
public:

  typedef typename Fb::Face_handle    Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
    typedef Delaunay_remove_tds_face_3_2<I, Fb2>           Other;
  };

  void set_info(const I &i) { inf = i; }

  const I& info() const { return inf; }

  // Handling the doubly connected list of faces
  // These functions should not be public, but private
  // and the tds should be declared friend.

  // Returns the sucessor
  Face_handle n() const {return _n;}

  // Returns the predecessor
  Face_handle p() const {return _p;}

  void set_p(Face_handle f) {_p = f;}

  void set_n(Face_handle f) {_n = f;}

private:
  // Dirty, so we should avoid it in another way
  Face_handle handle() const {
    Face_handle n = this->neighbor(0);
    int i = Triangulation_utils_3::ccw(n->index(this->vertex(1)));
    return n->neighbor(i);
  }

protected:
  // Remove this face from the list
  void remove_from_list() {
    // Q: Can't we be sure that there is always a predecessor
    // and successor?? This might pose a problem when we
    // remove the final tetrahedron, that is we have to
    // check whether that one still performs the
    // Surface::remove_vertex_3() method
    if(&*_p)
      _p->set_n(_n);

    if(&*_n)
      _n->set_p(_p);
  }

public:
  // Remove neighbor cw(i) and ccw(i) from the list
  void remove_neighbors_from_list(int i) {
    Face_handle n = this->neighbor(Triangulation_utils_3::cw(i));
    n->remove_from_list();
    n = this->neighbor(Triangulation_utils_3::ccw(i));
    n->remove_from_list();
  }

  // Marks edge i, that is marks one of the two half-edges,
  // namely the one in the face with the smaller address,
  // of the faces this and neighbor(i)
  //
  // Additionally, this face is then moved right behind face h,
  // because it is a candidate for an ear.
  void mark_edge(int i, Face_handle h) {
    Face_handle n = this->neighbor(i);
    if (n < handle()) {
      n->mark_halfedge(n->index(handle()));
      unmark_halfedge(i);
      h->move_after_this(n);
    } else {
      n->unmark_halfedge(n->index(handle()));
      mark_halfedge(i);
      if (h != handle())
	h->move_after_this(handle());
    }
  }

  // unmarks the two halfedges
  void unmark_edge(int i) {
    Face_handle n = this->neighbor(i);
    unmark_halfedge(i);
    int fi = n->index(handle());
    n->unmark_halfedge(fi);
  }

  // marks all edges adjacent to the face
  void mark_adjacent_edges() {
    for(int i = 0; i < 3; i++)
      mark_edge(i, handle());
  }

  bool is_halfedge_marked(int i) const {
    return _edge[i];
  }

  void set_edge(int i, bool b) {
    _edge[i] = b;
  }

protected:
  // Move face f after this.
  void move_after_this(Face_handle f) {
    if (_n == f)
      return;

    Face_handle p = f->p();
    Face_handle n = f->n();
    p->set_n(n);
    n->set_p(p);

    n = _n;
    _n = f;
    f->set_n(n);
    n->set_p(f);
    f->set_p(handle());
  }

  void unmark_halfedge(int i) {
    _edge[i] = false;
  }

  void mark_halfedge(int i) {
    _edge[i] = true;
  }

  I inf;
  Face_handle _p, _n;
  bool _edge[3];
};


// This class is used to represent the boundary of a hole in a polyhedron.
// It only implements a constructor, the rest is inherited

template <class T>
class Delaunay_remove_tds_3_2
  : public Triangulation_data_structure_2<
           Delaunay_remove_tds_vertex_3_2<typename T::Vertex_handle>,
	   Delaunay_remove_tds_face_3_2<typename T::Facet> >
{
public:
  typedef typename T::Facet Facet;
  typedef typename T::Cell_handle Cell_handle;
  typedef typename T::Vertex_handle Vertex_handle;

private:
  typedef Triangulation_data_structure_2<
            Delaunay_remove_tds_vertex_3_2<typename T::Vertex_handle>,
            Delaunay_remove_tds_face_3_2<typename T::Facet> > TDS_2;

public:
  typedef typename TDS_2::Face_iterator Face_iterator;
  typedef typename TDS_2::Face          Face_3_2;
  typedef typename TDS_2::Vertex_handle Vertex_handle_3_2;
  typedef typename TDS_2::Face_handle   Face_handle_3_2;

  using TDS_2::create_face;
  using TDS_2::create_vertex;
  using TDS_2::cw;
  using TDS_2::ccw;
  using TDS_2::set_dimension;
  using TDS_2::faces_begin;
  using TDS_2::faces_end;

private:

  struct Halfedge
    : public std::pair<Vertex_handle_3_2, Vertex_handle_3_2>
  {
    typedef std::pair<Vertex_handle_3_2, Vertex_handle_3_2>   Pair;

    Face_handle_3_2    third;
    int                fourth;

    Halfedge(Vertex_handle_3_2 a, Vertex_handle_3_2 b,
	     Face_handle_3_2 c, int d)
      : Pair(a, b), third(c), fourth(d) {}
  };

public:

  Delaunay_remove_tds_3_2(const std::vector<Facet> & boundhole );

  void remove_degree_3(Vertex_handle_3_2 v, Face_handle_3_2 f)
  {
    int i = f->index(v);
    // As f->neighbor(cw(i)) and f->neighbor(ccw(i)) will be removed,
    // we first remove it from the list we maintain.
    f->remove_neighbors_from_list(i);
    // we are ready to call the method of the base class
    TDS_2::remove_degree_3(v,f);
  }
};

template <class T>
Delaunay_remove_tds_3_2<T>::
Delaunay_remove_tds_3_2(const std::vector<Facet> & boundhole)
{
// FIXME : similar to operator>>(), isn't it ?  Should we try to factorize ?
    std::vector<Halfedge> halfedges;
    halfedges.reserve(3*boundhole.size());

    std::map<Vertex_handle, Vertex_handle_3_2>  vertex_map;
    typename std::map<Vertex_handle, Vertex_handle_3_2>::iterator map_it;

    for(typename std::vector<Facet>::const_iterator fit = boundhole.begin();
	    fit != boundhole.end(); ++fit)
    {
      Face_handle_3_2 f = create_face();

      Facet facet = *fit;
      Cell_handle h = facet.first;
      int k = facet.second;

      // All 2d faces must have the same orientation,
      // so we need a mapping from 3d to 2d indices.
      // Furthermore the triangles are seen "from the other side"
      int i0 = 0, i1 = 2, i2 = 3;
      switch (k) {
      case 0:  i0 = 1; break;
      case 1:  i1 = 3; i2 = 2; break;
      case 2:  i1 = 1; break;
      default: i2 = 1;
      }

      // We create as many 2d vertices as there are 3d vertices.
      Vertex_handle_3_2 v0, v1, v2;

      Vertex_handle w0 = h->vertex(i0);

      if((map_it = vertex_map.find(w0)) == vertex_map.end()) {
	v0 = create_vertex();
	v0->set_info(w0);
	vertex_map.insert(std::make_pair(w0, v0));
      } else {
	v0 = map_it->second;
      }

      Vertex_handle w1 = h->vertex(i1);
      if((map_it = vertex_map.find(w1)) == vertex_map.end()) {
	v1 = create_vertex();
	v1->set_info(w1);
	vertex_map.insert(std::make_pair(w1, v1));
      } else {
	v1 = map_it->second;
      }

      Vertex_handle w2 = h->vertex(i2);
      if((map_it = vertex_map.find(w2)) == vertex_map.end()) {
	v2 = create_vertex();
	v2->set_info(w2);
	vertex_map.insert(std::make_pair(w2, v2));
      } else {
	v2 = map_it->second;
      }

      v0->set_face(f);
      v1->set_face(f);
      v2->set_face(f);

      f->set_vertices(v0, v1, v2);
      f->set_info(facet);

      for(int j = 0; j < 3; j++) {
	Vertex_handle_3_2 v = f->vertex(j);
	Vertex_handle_3_2 w = f->vertex(cw(j));
	halfedges.push_back((v < w) ? Halfedge(v, w, f, ccw(j))
	                            : Halfedge(w, v, f, ccw(j)));
      }
    }

    std::sort(halfedges.begin(), halfedges.end());

    // The halfedges that are opposite to each other are neighbors
    // in the sorted list.

    for(typename std::vector<Halfedge>::iterator it = halfedges.begin();
	it != halfedges.end(); ++it) {
      Halfedge e1 = *it;
      ++it;
      Halfedge e2 = *it;
      e1.third->set_neighbor(e1.fourth, e2.third);
      e2.third->set_neighbor(e2.fourth, e1.third);
    }
    // The TDS cannot know that it is 2D because we constructed it
    // with advanced functions
    set_dimension(2);

    Face_iterator fit2 = faces_begin();
    Face_handle_3_2 first = fit2, f = fit2;
    for(int i = 0; i < 3; i++) {
	// we mark an edge only on one side
	f->set_edge(i, f < f->neighbor(i));
      }
    for(++fit2; fit2 != faces_end(); ++fit2) {
      f->set_n(fit2);
      fit2->set_p(f);
      f = fit2;
      for(int i = 0; i < 3; i++) {
	// we mark an edge only on one side
	f->set_edge(i, f < f->neighbor(i));
      }
    }
    // f points to the last face
    f->set_n(first);
    first->set_p(f);
}

}} // namespace CGAL::internal

#endif // CGAL_INTERNAL_DELAUNAY_REMOVE_TDS_3_H
