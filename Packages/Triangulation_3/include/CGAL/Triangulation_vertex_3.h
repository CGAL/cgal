// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Triangulation_vertex_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_VERTEX_3_H
#define CGAL_TRIANGULATION_VERTEX_3_H

#include <CGAL/Pointer.h>
//#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_circulators_3.h>
#include <CGAL/Triangulation_cell_3.h>
#include <CGAL/Triangulation_vertex_3.h>
#include <CGAL/Triangulation_handles_3.h>

#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds >
class Triangulation_cell_3;

template < class Gt, class Tds >
class Triangulation_vertex_handle_3;

template < class Gt, class Tds >
class Triangulation_cell_handle_3;

template<class Gt, class Tds >
class Triangulation_vertex_3
  : public Tds::Vertex
{
public:
  
  typedef typename Gt::Point_3 Point;

  typedef typename Tds::Vertex Vtds;
  typedef typename Tds::Cell Ctds;

  typedef Triangulation_cell_3<Gt,Tds> Cell;
  
  typedef Triangulation_vertex_handle_3<Gt,Tds> Vertex_handle;
  typedef Triangulation_cell_handle_3<Gt,Tds> Cell_handle;
  
  Triangulation_vertex_3()
     : Vtds()
  {}

  Triangulation_vertex_3(const Point & p)
    :  Vtds(p)
  {}
    
  Triangulation_vertex_3(const Point & p, Cell_handle c)
    :  Vtds(p, &(*c))
  {}

  Triangulation_vertex_3(Cell_handle c)
    :  Vtds(&(*c))
  {}

  void set_cell(Cell_handle c)
  {
    Vtds::set_cell(&(*c));
  }
    
  void set_point(const Point & p)
  {
    Vtds::set_point(p);
  }
    
  Cell_handle cell() const
  {
    return (Cell *) Vtds::cell();
  }
        
  Vertex_handle handle()
  {
    return Vertex_handle(this);
  }

  bool is_valid(bool verbose = false, int level = 0) const
  {
    return Vtds::is_valid(verbose,level);
  }
//   Vertex_circulator incident_vertices()
//   {
//     return Vertex_circulator(handle(), cell());
//   }

//   Vertex_circulator incident_vertices(const Cell_handle& c)
//   {
//     return Vertex_circulator(handle(), c);
//   } 

//   Cell_circulator incident_cells()
//   {
//     return Cell_circulator(handle(), cell());
//   }
    
//   Cell_circulator incident_cells(const Cell_handle& c)
//   {
//     return Cell_circulator(handle(), c);
//   }
    
//   Edge_circulator incident_edges()
//   {
//     return Edge_circulator(handle(), cell());
//   }
 
//   Edge_circulator incident_edges(const Cell_handle& c)
//   {
//     return Edge_circulator(handle(), c);
//   }
    
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_VERTEX_3_H
