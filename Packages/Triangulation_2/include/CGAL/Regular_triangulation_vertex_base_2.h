// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Regular_triangulation_vertex_base_2.h
// package       : Triangulation_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_2_H
#define CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_2_H

#include <CGAL/Triangulation_vertex_base_2.h>

CGAL_BEGIN_NAMESPACE

template < class GT,
           class Vbb = Triangulation_vertex_base_2<GT> >
class Regular_triangulation_vertex_base_2 
  :   public Vbb
{
  typedef typename Vbb::Triangulation_data_structure     TDS;
  typedef Vbb                                            Base; 
public:
  typedef typename Base::Point                Point;
  typedef TDS                                 Triangulation_data_structure;
  typedef typename TDS::Face_handle           Face_handle;
  typedef typename TDS::Vertex_handle         Vertex_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vbb::template Rebind_TDS<TDS2>::Other    Vb2;
    typedef Regular_triangulation_vertex_base_2<GT, Vb2>      Other;
  };

  Regular_triangulation_vertex_base_2 ()
    : Base(), _hidden(false)     {}
    
  Regular_triangulation_vertex_base_2(const Point & p, Face_handle f = NULL)
    :  Base(p, f), _hidden(false)
    {}

  void set_hidden(bool b) { _hidden = b; }
  bool is_hidden() { return _hidden ;}
 
private:
  bool _hidden;

};

CGAL_END_NAMESPACE

#endif //CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_2_H
