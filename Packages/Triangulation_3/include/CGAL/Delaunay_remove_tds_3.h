// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Delaunay_remove_tds_3.h
// revision      : $Revision$
//
// author(s)     : Monique Teillaud, Andreas Fabri
//
// coordinator   : INRIA Sophia Antipolis
//
// ============================================================================

#ifndef CGAL_DELAUNAY_REMOVE_TDS_3_H
#define CGAL_DELAUNAY_REMOVE_TDS_3_H

#include <CGAL/basic.h>

#include <map>
#include <CGAL/quadruple.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_using_list_2.h>

CGAL_BEGIN_NAMESPACE

// Instead of deriving from Triangulation_vertex_base_2<Gt> we copy the code
// Then we don't have the Point_2 problem

template <class Gt, class I>
class Delaunay_remove_tds_vertex_3_2  {
public :
  typedef typename Gt::Point_3  Point;

  Delaunay_remove_tds_vertex_3_2() 
    : _f(NULL) 
    {}

  Delaunay_remove_tds_vertex_3_2(const I & i)
    : _info(i)
    {}

  Delaunay_remove_tds_vertex_3_2(const Point & p, void * f = NULL) 
    : _f(f) 
    {}
  
  void* 
  face() const {     return _f;
  }

  void 
  set_face(void* f) {
    _f = f ;
  }

  bool
  is_valid(bool /* verbose */ = false, int /* level */ = 0) const {
      return true;
    }

  void 
  set_info(I i) { 
    _info = i;
  }

  I 
  info(){ 
    return _info; 
  }
  
  Point 
  point() {
    return _info->point();
  }

private:
  I _info;
  void * _f;
};


/* We derive the face class, 
   - because we want an additional pointer,
   - because we want to manage the order of the list of the faces
     in the triangulation datastructure,
   - because we want to look at each 'edge' only once. If two 
     triangles are neighbors, then the triangle with the smaller 
     address has the edge.
*/


template < class Gt,class I >
class Delaunay_remove_tds_face_3_2 :public Triangulation_face_base_2<Gt> {

public:
  Delaunay_remove_tds_face_3_2() {}

  Delaunay_remove_tds_face_3_2(void* v0, void* v1, void* v2) :
    Triangulation_face_base_2<void*>(v0, v1, v2) {}

  Delaunay_remove_tds_face_3_2(void* v0, void* v1, void* v2, 
			       void* n0, void* n1, void* n2) :
    Triangulation_face_base_2<void*>(v0, v1, v2, n0, n1, n2) {}

  void 
  set_info(I i) { 
    inf = i;    
  }

  I info(){ 
    return inf; 
  }

  //Handling the doubly connected list of faces
  void* 
  p() const {return _p;}

  void* 
  n() const {return _n;}

  void  
  set_p(void* f) {_p = f;};

  void
  set_n(void* f) {_n = f;};

  // Remove this face from the list
  void 
  unlink() {
    if(_p) {
      ((Delaunay_remove_tds_face_3_2*)_p)->set_n(_n);
    }
    if(_n) {
      ((Delaunay_remove_tds_face_3_2*)_n)->set_p(_p);
    }
  }

  // Remove neighbor i from the list
  void 
  unlink(int i) {
    Delaunay_remove_tds_face_3_2 * n = 
      (Delaunay_remove_tds_face_3_2*)neighbor(cw(i));
    n->unlink();
    n = (Delaunay_remove_tds_face_3_2*)neighbor(ccw(i));
    n->unlink();
  }

  // Move a given face after this
private:
  void 
  move_after_this(Delaunay_remove_tds_face_3_2* f) {
    if(_n == f) {
      return;
    }
    Delaunay_remove_tds_face_3_2 *p = (Delaunay_remove_tds_face_3_2*)f->p();
    Delaunay_remove_tds_face_3_2 *n = (Delaunay_remove_tds_face_3_2*)f->n();
    p->set_n(n);
    n->set_p(p);
    
    n = (Delaunay_remove_tds_face_3_2*)_n;
    _n = f;
    f->set_n(n);
    n->set_p(f);
    f->set_p(this);
  }

public:
  // Note that when we set an edge, the face gets moved
  // so that it gets considered later.
  void
  set_edge(int i, Delaunay_remove_tds_face_3_2* f) {
    Delaunay_remove_tds_face_3_2 *n = 
      (Delaunay_remove_tds_face_3_2*)neighbor(i);
    if(n < this) {
      n->set_edge(n->face_index(this), true);
      set_edge(i, false);
      f->move_after_this(n);
    } else {
      n->set_edge(n->face_index(this), false);
      set_edge(i, true);
    }      
  }

  void 
  set_edge(int i) {
    set_edge(i, this);
  }

  void 
  set_edge() {
    for(int i = 0; i < 3; i++) {
      set_edge(i, this);
    }
  }

  void 
  set_edge(int i, bool b) {
    _edge[i] = b;
  }


  bool 
  edge(int i) const {
    return _edge[i];
  }
    
private:
  I inf;
 
  void* _p;
  void* _n;

  bool _edge[3];

};


// The compare function is needed to find opposite edges. In fact it sorts
// a tuple by the first and second component. Maybe this is of general interest.
template < class H>
class Delaunay_remove_tds_halfedge_compare_3_2 {
public:
  bool operator()(const H & x, const H & y) {
    return ( x.first < y.first || 
	     ( (x.first == y.first) && (x.second < y.second)) );
  }
};


// This class is used to represent the boundary of a hole in a polyhedron.
// It only implements a constructor, the rest is inherited

template <class T>
class Delaunay_remove_tds_3_2 
: 
  public Triangulation_data_structure_using_list_2< 
           Delaunay_remove_tds_vertex_3_2<typename T::Geom_traits,
                                          typename T::Vertex*>,
	   Delaunay_remove_tds_face_3_2<typename T::Geom_traits, 
	                                typename T::Facet> > 
{
public:
  typedef typename T::Facet Facet;
  typedef typename T::Cell_handle Cell_handle;
  typedef typename T::Vertex_handle Vertex_handle;
  typedef typename T::Vertex Vertex;

private:
  typedef Triangulation_data_structure_using_list_2< 
            Delaunay_remove_tds_vertex_3_2<typename T::Geom_traits,
	                                   typename T::Vertex*>,
            Delaunay_remove_tds_face_3_2<typename T::Geom_traits,
	                                 typename T::Facet> > TDSUL2;

public:
  typedef typename TDSUL2::Vertex Vertex_3_2;
  typedef typename TDSUL2::Face Face_3_2;
  typedef typename TDSUL2::Face_iterator Face_iterator;
  typedef quadruple<void*, void*, Face_3_2*, int> Halfedge;



  Delaunay_remove_tds_3_2( std::list<Facet> & boundhole ) {

    typedef typename std::list<Facet>::iterator Facet_iterator;

    int size = boundhole.size();
 
    std::vector<Halfedge> halfedges(3*size);
    int i = 0;
    std::map<Vertex*, Vertex_3_2*>  vertex_map;
    typename std::map<Vertex*, Vertex_3_2*>::iterator map_it;

    for(Facet_iterator fit = boundhole.begin() ; 
	fit != boundhole.end(); 
	++fit) {
      Face_3_2 * f = create_face();

      Facet facet = *fit;
      Cell_handle h = facet.first;
      int k = facet.second;
    
      // All 2d faces must have the same orientation, 
      // so we need a mapping from 3d to 2d indices.
      // Furthermore the triangles are seen "from the other side"
      int i0 = 0, i1 = 2, i2 = 3;
      switch (k) {
      case 0: i0 = 1;   break;
      case 1:  i1 = 3; i2 = 2; break;
      case 2:   i1 = 1; break;
      default:  i2 = 1;
      }

      // We create as many 2d vertices as there are 3d vertices.
      Vertex_3_2 *v0, *v1, *v2;

      Vertex *w0 = handle2pointer(h->vertex(i0));

      if((map_it = vertex_map.find(w0)) == vertex_map.end()) {
	v0 = create_vertex();
	v0->set_info(w0);
	vertex_map.insert(std::make_pair(w0, v0));
      } else {
	v0 = map_it->second;
      }
      
      Vertex *w1 = handle2pointer(h->vertex(i1));
      if((map_it = vertex_map.find(w1)) == vertex_map.end()) {
	v1 = create_vertex();
	v1->set_info(w1);
	vertex_map.insert(std::make_pair(w1, v1));
      } else {
	v1 = map_it->second;
      }
      
      Vertex *w2 = handle2pointer(h->vertex(i2));
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
      ///

      for(int j = 0; j < 3; j++) {
	void* v = f->vertex(j); 
	void* w = f->vertex(cw(j));
	halfedges[i] = (v < w)? make_quadruple(v, w, f, ccw(j))
	                      : make_quadruple(w, v, f, ccw(j));
	i++;
      }
    }

    std::sort(halfedges.begin(), 
	      halfedges.end(), 
	      Delaunay_remove_tds_halfedge_compare_3_2<Halfedge>());


    // The halfedges that are oppsoite to each other are neighbor
    // in the sorted list. 

    for(typename std::vector<Halfedge>::iterator it = halfedges.begin();
	it != halfedges.end();
	++it) {
      Halfedge e1 = *it;
      ++it;
      Halfedge e2 = *it;
      e1.third->set_neighbor(e1.fourth, e2.third);
      e2.third->set_neighbor(e2.fourth, e1.third);
    }
    // The TDS cannot know that it is 2D because we constructed it 
    // with advanced functions
    set_dimension(2);

    Face_3_2 * dummy = new Face_3_2();;

    Face_3_2 *f = dummy;

    for( Face_iterator fit = faces_begin(); fit != faces_end(); ++fit) {
      f->set_n(&(*fit));
      fit->set_p(f);
      f = &(*fit);
      for(int i = 0; i < 3; i++) {
	f->set_edge(i, (f < (f->neighbor(i))));
      }
    }
    // f points to the last face
    f->set_n(dummy->n());
    ((Face_3_2*)dummy->n())->set_p(f);
    delete(dummy);
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_DELAUNAY_REMOVE_TDS_3_H
