// Copyright (c)  2005, 2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Raphaelle Chaine  (Raphaelle.Chaine@sophia.inria.fr, raphaelle.chaine@liris.cnrs.fr)
//                 Andreas Fabri (GeometryFactory)                 
//

#ifndef CGAL_TFS_VERTEX_BASE_3_H
#define CGAL_TFS_VERTEX_BASE_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

template < class Vbb >
class TFS_vertex_base_3
  : public Vbb
{
  typedef Vbb                                           Base;
  typedef typename Base::Triangulation_data_structure   Tds;
public:
  typedef typename Tds::Vertex_handle                   Vertex_handle;
  typedef typename Tds::Cell_handle                     Cell_handle;
  typedef typename Base::Point                          Point;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vbb::template Rebind_TDS<TDS2>::Other      Vb2;
    typedef TFS_vertex_base_3<Vb2>          Other;
  };

  TFS_vertex_base_3()
    : Base(), slice_index(-1) {}

  TFS_vertex_base_3(const Point & p, Cell_handle f,int num=-1)
    : Base(p,f), slice_index(num) {}

  TFS_vertex_base_3(const Point & p,int num=-1)
    : Base(p), slice_index(num) {}

  int slice() const { return slice_index; }
  void set_slice_index(int num)   { slice_index=num; }

private:
  int slice_index;
};

CGAL_END_NAMESPACE

#endif // CGAL_TFS_VERTEX_BASE_3_H
