// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : include/CGAL/Apollonius_graph_hierarchy_vertex_base_2.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================


#ifndef CGAL_APOLLONIUS_GRAPH_HIERARCHY_VERTEX_BASE_2_H
#define CGAL_APOLLONIUS_GRAPH_HIERARCHY_VERTEX_BASE_2_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template < class Vbb >
class Apollonius_graph_hierarchy_vertex_base_2
 : public Vbb
{
  typedef Vbb                                              Base;
  typedef typename Base::Apollonius_graph_data_structure   Agds;

public:
  typedef typename Base::Site_2             Site_2;
  typedef Agds                              Apollonius_graph_data_structure;
  typedef typename Agds::Vertex_handle      Vertex_handle;
  typedef typename Agds::Face_handle        Face_handle;

  template < typename AGDS2 >
  struct Rebind_TDS {
    typedef typename Vbb::template Rebind_TDS<AGDS2>::Other   Vb2;
    typedef Apollonius_graph_hierarchy_vertex_base_2<Vb2>     Other;
  };

  Apollonius_graph_hierarchy_vertex_base_2()
    : Base(), _up(NULL), _down(NULL)
    {}
  Apollonius_graph_hierarchy_vertex_base_2(const Site_2& p,
					   Face_handle f)
    : Base(p,f), _up(NULL), _down(NULL)
    {}
  Apollonius_graph_hierarchy_vertex_base_2(const Site_2& p)
    : Base(p), _up(NULL), _down(NULL)
    {}

  Vertex_handle up() {return _up;}
  Vertex_handle down() {return _down;}
  void set_up(Vertex_handle u) {_up=u;}
  void set_down(Vertex_handle d) {if (this) _down=d;}


 private:
  Vertex_handle  _up;    // same vertex one level above
  Vertex_handle  _down;  // same vertex one level below
};

CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_HIERARCHY_VERTEX_BASE_2_H
