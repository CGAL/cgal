
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
// file          : include/CGAL/Triangulation_circulators_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (Mariette Yvinec)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_CIRCULATORS_3_H
#define CGAL_TRIANGULATION_CIRCULATORS_3_H

#include <pair.h>
#include <CGAL/triple.h>
#include <CGAL/circulator.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Triangulation_ds_circulators_3.h>

template < class Gt, class Tds >
class CGAL_Triangulation_3;

template < class Gt, class Tds >
class CGAL_Triangulation_cell_3;

template < class Gt, class Tds>
class CGAL_Triangulation_cell_circulator_3
  : public CGAL_Bidirectional_circulator_base<CGAL_Triangulation_cell_3<Gt,Tds>, ptrdiff_t, size_t>
{
public:
  typedef typename Tds::Cell Ctds;
  typedef typename Tds::Cell_circulator Circulator_base;

  typedef CGAL_Triangulation_cell_3<Gt,Tds> Cell;
  typedef CGAL_Triangulation_vertex_3<Gt,Tds> Vertex;
  typedef typename Vertex::Vertex_handle Vertex_handle;
  typedef typename Cell::Cell_handle Cell_handle;
  typedef CGAL_Triangulation_3<Gt,Tds> Triangulation_3;
  typedef typename Triangulation_3::Edge Edge;

  typedef CGAL_Triangulation_cell_circulator_3<Gt,Tds> Cell_circulator;

  CGAL_Triangulation_cell_circulator_3()
    : _cb(), _tr(NULL)
    {}

  CGAL_Triangulation_cell_circulator_3(Triangulation_3 * tr, Edge e)
    : _cb( &(tr->_tds), CGAL_make_triple( (Ctds *) &(*(e.first)), e.second, e.third ) ), _tr(tr)
    {}

   CGAL_Triangulation_cell_circulator_3(Triangulation_3 * tr, Edge e, Cell_handle c)
    : _cb( &(tr->_tds), 
	   CGAL_make_triple( (Ctds *) &(*(e.first)), e.second, e.third ),
	   (Ctds *) &(*c) ), 
      _tr(tr)
    {}
 
  CGAL_Triangulation_cell_circulator_3(const Cell_circulator & ccir)
    : _cb(ccir._cb), _tr(ccir._tr)
    {}

  Cell_circulator&
  operator=(const Cell_circulator & ccir)
    {
      _cb = ccir._cb;
      _tr = ccir._tr;
      return *this;
    }
  
  bool
  operator==(const Cell_circulator & ccir) const
  {
    return ( _cb == ccir._cb);
  }

  bool
  operator!=(const Cell_circulator & ccir)
  {
    return ( !(*this == ccir) );
  }

  Cell_circulator &
  operator++()
  {
    ++_cb;
    return *this;
  }

  Cell_circulator &
  operator--()
  {
    --_cb;
    return *this;
  }

  Cell_circulator
  operator++(int)
  {
    Cell_circulator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Cell_circulator
  operator--(int)
  {
    Cell_circulator tmp(*this);
    --(*this);
    return tmp;
  }

    inline Cell & operator*() const
  {
    return (Cell &)(*_cb);
  }

  inline Cell* operator->() const
  {
    return (Cell*)( &(*_cb) );
  }

private: 
  Circulator_base _cb;
  Triangulation_3 * _tr;
};

#endif
