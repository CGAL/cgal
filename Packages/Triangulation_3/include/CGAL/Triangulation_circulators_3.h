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
// file          : include/CGAL/Triangulation_circulators_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_CIRCULATORS_3_H
#define CGAL_TRIANGULATION_CIRCULATORS_3_H

#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/circulator.h>
#include <CGAL/Triangulation_ds_circulators_3.h>
#include <CGAL/Triangulation_vertex_3.h>
#include <CGAL/Triangulation_cell_3.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds >
class Triangulation_3;

template < class Gt, class Tds>
class Triangulation_cell_circulator_3
  : public Bidirectional_circulator_base<Triangulation_cell_3<Gt,Tds>, 
    ptrdiff_t, size_t>
{
public:
  typedef typename Tds::Cell Ctds;
  typedef typename Tds::Cell_circulator Circulator_base;

  typedef Triangulation_3<Gt,Tds> Triangulation;

  typedef typename Triangulation::Cell Cell;
  typedef typename Triangulation::Vertex Vertex;
  typedef typename Triangulation::Edge Edge;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Cell_handle Cell_handle;

  typedef Triangulation_cell_circulator_3<Gt,Tds> Cell_circulator;

  Triangulation_cell_circulator_3()
    : _cb(), _tr(NULL)
    {}

  Triangulation_cell_circulator_3(const Triangulation * tr, Cell_handle c,
				  int s, int t)
    : _cb( &(const_cast<Triangulation *>(tr)->_tds), (Ctds *) &(*c), s, t ),
      _tr(const_cast<Triangulation *>(tr))
    {}

  Triangulation_cell_circulator_3(const Triangulation * tr, const Edge & e)
    : _cb( &(const_cast<Triangulation *>(tr)->_tds), (Ctds *) &(*(e.first)),
	    e.second, e.third ), _tr(const_cast<Triangulation *>(tr))
    {}

   Triangulation_cell_circulator_3(const Triangulation * tr, 
				   Cell_handle c, int s, int t,
				   Cell_handle start)
    : _cb( &(const_cast<Triangulation *>(tr)->_tds), 
	   (Ctds *) &(*c), s, t, (Ctds *) &(*start) ),
      _tr(const_cast<Triangulation *>(tr))
    {}

   Triangulation_cell_circulator_3(const Triangulation * tr, const Edge & e, 
				   Cell_handle start)
    : _cb( &(const_cast<Triangulation *>(tr)->_tds), 
	   (Ctds *) &(*(e.first)), e.second, e.third, 
	   (Ctds *) &(*start) ),
      _tr(const_cast<Triangulation *>(tr))
    {}
 
  Triangulation_cell_circulator_3(const Cell_circulator & ccir)
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

  Cell & operator*() const
  {
    return (Cell &)(*_cb);
  }

  Cell* operator->() const
  {
    return (Cell*)( &(*_cb) );
  }

private: 
  Circulator_base _cb;
  Triangulation * _tr;
};

template < class Gt, class Tds>
class Triangulation_facet_circulator_3
  : public Bidirectional_circulator_base<typename 
                                         Triangulation_3<Gt,Tds>::Facet,
    ptrdiff_t, size_t>
{
public:
  typedef typename Tds::Cell Ctds;
  typedef typename Tds::Facet_circulator Circulator_base;

  typedef Triangulation_3<Gt,Tds> Triangulation;

  typedef typename Triangulation::Cell Cell;
  typedef typename Triangulation::Vertex Vertex;
  typedef typename Triangulation::Edge Edge;
  typedef typename Triangulation::Facet Facet;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Cell_handle Cell_handle;

  typedef Triangulation_facet_circulator_3<Gt,Tds> Facet_circulator;

  Triangulation_facet_circulator_3()
    : _cb(), _tr(NULL)
    {}

  Triangulation_facet_circulator_3(const Triangulation * tr, Cell_handle c,
				  int s, int t)
    : _cb( &(const_cast<Triangulation *>(tr)->_tds), (Ctds *) &(*c), s, t ),
      _tr(const_cast<Triangulation *>(tr))
    {}

  Triangulation_facet_circulator_3(const Triangulation * tr, const Edge & e)
    : _cb( &(const_cast<Triangulation *>(tr)->_tds), (Ctds *) &(*(e.first)),
	    e.second, e.third ), _tr(const_cast<Triangulation *>(tr))
    {}

   Triangulation_facet_circulator_3(const Triangulation * tr, 
				   Cell_handle c, int s, int t,
				   const Facet & start)
    : _cb( &(const_cast<Triangulation *>(tr)->_tds),
	   (Ctds *) &(*c), s, t,
	   std::make_pair((Ctds *) &(*start.first), start.second) ),
      _tr(const_cast<Triangulation *>(tr))
    {}

   Triangulation_facet_circulator_3(const Triangulation * tr, const Edge & e, 
				   const Facet & start)
    : _cb( &(const_cast<Triangulation *>(tr)->_tds),
	   (Ctds *) &(*(e.first)), e.second, e.third, 
	   std::make_pair((Ctds *) &(*start.first), start.second) ),
      _tr(const_cast<Triangulation *>(tr))
    {}
 
   Triangulation_facet_circulator_3(const Triangulation * tr, 
				   Cell_handle c, int s, int t,
				   Cell_handle start, int f)
    : _cb( &(const_cast<Triangulation *>(tr)->_tds),
	   (Ctds *) &(*c), s, t,
	   (Ctds *) &(*start), f ),
      _tr(const_cast<Triangulation *>(tr))
    {}

   Triangulation_facet_circulator_3(const Triangulation * tr, const Edge & e, 
				   Cell_handle start, int f)
    : _cb( &(const_cast<Triangulation *>(tr)->_tds),
	   (Ctds *) &(*(e.first)), e.second, e.third, 
	   (Ctds *) &(*start), f ),
      _tr(const_cast<Triangulation *>(tr))
    {}
 
  Triangulation_facet_circulator_3(const Facet_circulator & ccir)
    : _cb(ccir._cb), _tr(ccir._tr)
    {}

  Facet_circulator&
  operator=(const Facet_circulator & ccir)
    {
      _cb = ccir._cb;
      _tr = ccir._tr;
      return *this;
    }
  
  bool
  operator==(const Facet_circulator & ccir) const
  {
    return ( _cb == ccir._cb);
  }

  bool
  operator!=(const Facet_circulator & ccir)
  {
    return ( !(*this == ccir) );
  }

  Facet_circulator &
  operator++()
  {
    ++_cb;
    return *this;
  }

  Facet_circulator &
  operator--()
  {
    --_cb;
    return *this;
  }

  Facet_circulator
  operator++(int)
  {
    Facet_circulator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Facet_circulator
  operator--(int)
  {
    Facet_circulator tmp(*this);
    --(*this);
    return tmp;
  }

  Facet operator*() const
  {
    return std::make_pair( Cell_handle( (Cell *) ((*_cb).first) ),
			   (*_cb).second  ) ;
  }

private: 
  Circulator_base _cb;
  Triangulation * _tr;
};

CGAL_END_NAMESPACE

#endif  // CGAL_TRIANGULATION_CIRCULATORS_3_H
