
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
// file          : include/CGAL/Triangulation_circulators.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_CIRCULATORS_H
#define CGAL_TRIANGULATION_CIRCULATORS_H

#include <pair.h>
#include <CGAL/triple.h>
#include <CGAL/circulator.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names.h>
#include <CGAL/Triangulation_ds_circulators.h>

template < class Gt, class Tds >
class CGAL_Triangulation;

template < class Gt, class Tds >
class CGAL_Triangulation_cell;

template < class Gt, class Tds>
class CGAL_Triangulation_cell_circulator
  : public CGAL_Bidirectional_circulator_base<CGAL_Triangulation_cell<Gt,Tds>, ptrdiff_t, size_t>
{
public:
  typedef typename Tds::Cell Ctds;
  typedef typename Tds::Cell_circulator Circulator_base;

  typedef CGAL_Triangulation_cell<Gt,Tds> Cell;
  typedef CGAL_Triangulation_vertex<Gt,Tds> Vertex;
  typedef typename Vertex::Vertex_handle Vertex_handle;
  typedef typename Cell::Cell_handle Cell_handle;
  typedef CGAL_Triangulation<Gt,Tds> Triangulation;
  typedef typename Triangulation::Edge Edge;

  typedef CGAL_Triangulation_cell_circulator<Gt,Tds> Cell_circulator;

  CGAL_Triangulation_cell_circulator()
    : _cb(), _tr()
    {}

  CGAL_Triangulation_cell_circulator(Triangulation * tr, Edge e)
    : _cb( &(tr->_tds), e ), _tr(tr)
    {}
  
  CGAL_Triangulation_cell_circulator(const Cell_circulator & ccir)
    : _cb(ccir._cb), _tr(cit._tr)
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
    return ( _cb == cit._cb);
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
  }

  Cell_circulator &
  operator--()
  {
    --_cb;
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
  Triangulation * _tr;
};

#endif
