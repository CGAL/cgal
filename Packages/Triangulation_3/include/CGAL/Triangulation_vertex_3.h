// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : include/CGAL/Triangulation_vertex.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_VERTEX_H
#define CGAL_TRIANGULATION_VERTEX_H

#include <CGAL/Pointer.h>
#include <CGAL/Triangulation_short_names.h>
#include <CGAL/Triangulation_data_structure.h>
#include <CGAL/Triangulation_circulators.h>


template < class Gt, class Tds >
class CGAL_Triangulation_cell;

template < class Gt, class Tds >
class CGAL_Triangulation_vertex_handle;

template < class Gt, class Tds >
class CGAL_Triangulation_cell_handle;

// template<class Gt,  class Tds>
// class CGAL_Triangulation_cell_circulator;

// template<class Gt,  class Tds>
// class CGAL_Triangulation_vertex_circulator;

// template<class Gt, class Tds>
// class CGAL_Triangulation_edge_circulator;

template<class Gt, class Tds >
class CGAL_Triangulation_vertex
  : public Tds::Vertex
{
public:
  
  typedef typename Gt::Point Point;

  typedef typename Tds::Vertex Vtds;
  typedef typename Tds::Cell Ctds;

  typedef CGAL_Triangulation_cell<Gt,Tds> Cell;
  typedef CGAL_Triangulation_vertex<Gt,Tds> Vertex;
  
  typedef CGAL_Triangulation_vertex_handle<Gt,Tds> Vertex_handle;
  typedef CGAL_Triangulation_cell_handle<Gt,Tds> Cell_handle;
  //  typedef pair<Cell_handle, int, int>     Edge;

  //  typedef CGAL_Triangulation_cell_circulator<Gt,Tds>      Cell_circulator;
  //  typedef CGAL_Triangulation_edge_circulator<Gt,Tds>      Edge_circulator;
  //  typedef CGAL_Triangulation_vertex_circulator<Gt,Tds>    Vertex_circulator;

  
  CGAL_Triangulation_vertex()
     : Vtds()
  {}

  CGAL_Triangulation_vertex(const Point & p)
    :  Vtds(p)
  {}
    
//   CGAL_Triangulation_vertex(const Point & p, Cell_handle c)
//     :  Vtds(p, &(*c))
//   {}

  inline void set_cell(const Cell_handle & c)
  {
    Vtds::set_cell(&(*c));
  }
    
  inline Cell_handle cell() const
  {
    return (Cell *) Vtds::cell();
  }
        
  inline Vertex_handle handle() const
  {
    return Vertex_handle(this);
  }

//   inline Vertex_circulator incident_vertices()
//   {
//     return Vertex_circulator(handle(), cell());
//   }

//   inline Vertex_circulator incident_vertices(const Cell_handle& c)
//   {
//     return Vertex_circulator(handle(), c);
//   } 

//   inline 
//   Cell_circulator incident_cells()
//   {
//     return Cell_circulator(handle(), cell());
//   }
    
//   inline Cell_circulator incident_cells(const Cell_handle& c)
//   {
//     return Cell_circulator(handle(), c);
//   }
    
//   inline Edge_circulator incident_edges()
//   {
//     return Edge_circulator(handle(), cell());
//   }
 
//   inline Edge_circulator incident_edges(const Cell_handle& c)
//   {
//     return Edge_circulator(handle(), c);
//   }
    
};
#endif CGAL_TRIANGULATION_H
