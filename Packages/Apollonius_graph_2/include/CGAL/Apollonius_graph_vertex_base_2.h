// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Apollonius_graph_vertex_base_2.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_APOLLONIUS_GRAPH_VERTEX_BASE_2_H
#define CGAL_APOLLONIUS_GRAPH_VERTEX_BASE_2_H

#include <list>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/triangulation_assertions.h>

CGAL_BEGIN_NAMESPACE


template <class Gt,
	  bool StoreHidden = true,
	  class Vb = Triangulation_ds_vertex_base_2<> >
class Apollonius_graph_vertex_base_2
  : public Vb
{
private:
  typedef typename Vb::Triangulation_data_structure   AGDS;
public:
  // TYPES
  //------
  typedef Gt                             Geom_traits;
  typedef Vb                             Base;
  typedef typename Gt::Weighted_point_2  Weighted_point_2;
  typedef AGDS	                         Apollonius_graph_data_structure;
  typedef typename AGDS::Face_handle     Face_handle;
  typedef typename AGDS::Vertex_handle   Vertex_handle;

  template < typename AGDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<AGDS2>::Other      Vb2;
    typedef Apollonius_graph_vertex_base_2<Gt,StoreHidden,Vb2>  Other;
  };


private:
  // local types
  typedef std::list<Weighted_point_2>         Container;

public:
  // TYPES (continued)
  //------------------
  typedef Container                        Hidden_weighted_point_container;
  typedef typename Container::iterator     Hidden_weighted_point_iterator;

public:
  // CREATION
  //---------
  Apollonius_graph_vertex_base_2() : Vb() {}
  Apollonius_graph_vertex_base_2(const Weighted_point_2& p) : Vb(), _p(p) {}
  Apollonius_graph_vertex_base_2(const Weighted_point_2& p,
				 Face_handle f)
    : Vb(f), _p(p) {}

  ~Apollonius_graph_vertex_base_2()
  {
    clear_hidden_weighted_point_container();
  }


  // ACCESS METHODS
  //---------------
  Weighted_point_2 point() const { return _p; }

  Face_handle face() const { return Vb::face(); }

  unsigned int number_of_hidden_weighted_points() const {
    return weighted_point_list.size();
  }

  Hidden_weighted_point_iterator hidden_weighted_points_begin() { 
    return weighted_point_list.begin();
  }

  Hidden_weighted_point_iterator hidden_weighted_points_end() {
    return weighted_point_list.end();
  }

public:
  // SETTING AND UNSETTING
  //----------------------
  void set_point(const Weighted_point_2& p) { _p = p; }


  void add_hidden_weighted_point(const Weighted_point_2& p)
  {
    if ( StoreHidden ) {
      weighted_point_list.push_back(p);
    }
  }

  void clear_hidden_weighted_point_container()
  {
    weighted_point_list.clear();
  }

public:
  // VALIDITY CHECK
  bool is_valid(bool verbose = false, int level = 0) const {
    return Vb::is_valid(verbose, level);
  }

private:
  // class variables
  Container weighted_point_list;
  Weighted_point_2 _p;
};

CGAL_END_NAMESPACE 

#endif // CGAL_APOLLONIUS_GRAPH_VERTEX_BASE_2_H
