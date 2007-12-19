// Copyright (c) 2005, 2006  INRIA Sophia-Antipolis (France).
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

/*! \file TFS_polyline_vertex_base_3.h
    \brief Vertex class for the triangulation from slices.
*/

#ifndef CGAL_POLYLINE_VERTEX_BASE_3_H
#define CGAL_POLYLINE_VERTEX_BASE_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

/*!
  This vertex class stores additional data needed by the reconstruction algorithm.
*/
template < class Vbb >
class TFS_polyline_vertex_base_3
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
    typedef TFS_polyline_vertex_base_3<Vb2>          Other;
  };


  TFS_polyline_vertex_base_3()
    : Base(), slice_index(-1), m_next(NULL)  {}

  TFS_polyline_vertex_base_3(const Point & p, Cell_handle f,int num=-1)
    : Base(p,f), slice_index(num), m_next(NULL) {}

  TFS_polyline_vertex_base_3(const Point & p,int num=-1)
    : Base(p), slice_index(num), m_next(NULL) {}

/*!
 * The slice of the vertex in the triangulation.
 */
  int slice() const 
  { return slice_index; }

  void set_slice(int num)   { slice_index=num; }

/**
 * The next vertex on a polygonal contour.
 */
  Vertex_handle next() {return m_next;}
  void set_next(Vertex_handle v) {m_next=v;}

private:
  int slice_index; //Number of the slice the vertex belongs to
  Vertex_handle m_next; // Next vertex in the polygonal contour

//OFF file production
private:
  int m_index;
public :

/**
 * Obtain the index of a vertex. It is in the user responsibility to initialize the index.
 */
  int index() const
   { return m_index;}

/**
 * Set the index of a vertex.
 */
  void set_index(int n)
    { m_index = n; }

};

CGAL_END_NAMESPACE

#endif // CGAL_POLYLINE_VERTEX_BASE_3_H
