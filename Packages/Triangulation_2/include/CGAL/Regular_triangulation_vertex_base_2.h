// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri


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
    
  Regular_triangulation_vertex_base_2(const Point & p) 
    : Base(p), _hidden(false)   {}

  Regular_triangulation_vertex_base_2(const Point & p, Face_handle)
    : Base(p, f), _hidden(false) {}

  void set_hidden(bool b) { _hidden = b; }
  bool is_hidden() { return _hidden ;}
 
private:
  bool _hidden;

};

CGAL_END_NAMESPACE

#endif //CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_2_H
