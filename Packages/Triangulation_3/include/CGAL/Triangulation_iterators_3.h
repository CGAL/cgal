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
// file          : include/CGAL/Triangulation_iterators_3.h
// revision      : $Revision$
//
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_ITERATORS_3_H
#define CGAL_TRIANGULATION_ITERATORS_3_H

#include <CGAL/basic.h>
#include <utility>
#include <iterator>

#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/triangulation_assertions.h>

CGAL_BEGIN_NAMESPACE

template < class Tr>
class Triangulation_cell_iterator_3
{
  typedef typename Tr::Triangulation_data_structure Tds;
  typedef typename Tds::Cell                        Ctds;
  typedef typename Tds::Cell_iterator               Iterator_base;

  typedef typename Tr::Cell            Cell;
  typedef typename Tr::Vertex          Vertex;
  typedef typename Tr::Cell_iterator   Cell_iterator;

public:
  typedef Cell                              value_type;
  typedef Cell *                            pointer;
  typedef Cell &                            reference;
  typedef std::size_t                       size_type;
  typedef std::ptrdiff_t                    difference_type;
  typedef std::bidirectional_iterator_tag   iterator_category;

  Triangulation_cell_iterator_3()
    : _ib(), _tr(NULL), _inf(true)
  {}
        
  Triangulation_cell_iterator_3(const Tr *tr, bool inf)
    : _ib( tr->tds().cells_begin()), _tr(tr), _inf(inf)
    { 
      if (! _inf) {
	while ( // ( _ib != _tr->tds().cells_end() ) &&
	       // useless : there must be at least one finite cell
	       // since precond in _ib : dimension == 3
	       _tr->is_infinite(&*_ib) )
	  { ++_ib; }
      }
    }
  
  // for past-end iterator
  // does not need to find a finite cell
  Triangulation_cell_iterator_3(const Tr *tr)
    : _ib( tr->tds().cells_end()), _tr(tr), _inf(true)
  {}

  bool
  operator==(const Cell_iterator & cit) const
  {
    if ( _tr != cit._tr ) 
      return false;

    if ( ( _ib == _tr->tds().cells_end() ) 
	 || ( cit._ib == _tr->tds().cells_end() ) ) 
      return ( _ib == cit._ib );

    return ( ( _ib == cit._ib ) && ( _inf == cit._inf ) );
  }

  bool
  operator!=(const Cell_iterator & cit) const
  {
    return ( !(*this == cit) );
  }

  Cell_iterator &
  operator++()
  {
    if (_inf) {
      ++_ib;
    }
    else {
      do {
	++_ib; 
      } while ( ( _ib != _tr->tds().cells_end() )
		&& _tr->is_infinite(&*_ib) );
    }
    return *this;   
  }

  Cell_iterator &
  operator--()
  {
    if (_inf) {
      --_ib;
    }
    else{
      do {
	--_ib;
      } while ( ( _ib != _tr->tds().cells_end() )
		&& _tr->is_infinite(&*_ib) );
    }
    return *this;   
  }

  Cell_iterator
  operator++(int)
  {
    Cell_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Cell_iterator
  operator--(int)
  {
    Cell_iterator tmp(*this);
    --(*this);
    return tmp;
  }
        
  Cell & operator*() const
  {
    return (Cell &)(*_ib);
  }

  Cell* operator->() const
  {
    return (Cell*)( &(*_ib) );
  }
     
private: 
  Iterator_base _ib;
  const Tr * _tr;
  bool _inf; // if _inf == true, traverses all cells
               // else only traverses finite cells
};


template <class Tr>
class Triangulation_vertex_iterator_3
{
  typedef typename Tr::Triangulation_data_structure Tds;
  typedef typename Tds::Cell                        Ctds;
  typedef typename Tds::Vertex_iterator             Iterator_base;

  typedef typename Tr::Cell            Cell;
  typedef typename Tr::Vertex          Vertex;
  typedef typename Tr::Vertex_iterator Vertex_iterator;  
  typedef typename Tr::Cell_handle     Cell_handle;
  typedef typename Tr::Cell_iterator   Cell_iterator;

public:
  typedef Vertex                value_type;
  typedef Vertex *              pointer;
  typedef Vertex &              reference;
  typedef std::size_t           size_type;
  typedef std::ptrdiff_t        difference_type;
  typedef std::bidirectional_iterator_tag   iterator_category;

  Triangulation_vertex_iterator_3()
    : _ib(), _tr(NULL), _inf(true)
  {}

  Triangulation_vertex_iterator_3(const Tr * tr, bool inf)
    : _ib(tr->tds().vertices_begin()), _tr(tr), _inf(inf)
    { 
      if (! _inf) {
	if ( _tr->is_infinite(&*_ib) ) {
	  ++_ib;
	}
      }
    }
        
  // for past-end iterator
  // does not need to find a finite cell
  Triangulation_vertex_iterator_3(const Tr * tr)
    : _ib(tr->tds().vertices_end()), _tr(tr)
  { }
       
  bool
  operator==(const Vertex_iterator & vi) const
  {
    if ( _tr != vi._tr ) 
      return false;

    if ( ( _ib == _tr->tds().vertices_end() ) 
	 || ( vi._ib == _tr->tds().vertices_end() ) ) 
      return ( _ib == vi._ib );

    return ( ( _ib == vi._ib ) && ( _inf == vi._inf ) );
  }

  bool
  operator!=(const Vertex_iterator & vi) const
  {
    return ( !(*this == vi) );
  }

  Vertex_iterator &
  operator++()
  {
    ++_ib;
    if (!_inf
	&& _ib != _tr->tds().vertices_end()
        && _tr->is_infinite(&*_ib) )
	++_ib;
    return *this;   
  }

  Vertex_iterator &
  operator--()
  {
    --_ib;
    if (!_inf && _tr->is_infinite(&*_ib) )
	--_ib;
    return *this;   
  }

  Vertex_iterator
  operator++(int)
  {
    Vertex_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Vertex_iterator
  operator--(int)
  {
    Vertex_iterator tmp(*this);
    --(*this);
    return tmp;
  }

  Vertex & operator*() const
  {
    return (Vertex &)(*_ib);
  }

  Vertex* operator->() const
  {
    return   (Vertex*)( &(*_ib) );
  }
     
private:
  Iterator_base   _ib;
  const Tr * _tr;
  bool _inf; // if _inf == true, traverses all vertices
               // else only traverses finite vertices
};


template <class Tr>
class Triangulation_edge_iterator_3
{
  typedef typename Tr::Triangulation_data_structure Tds;
  typedef typename Tr::Edge  Edge;
  typedef typename Tr::Cell  Cell;
  typedef typename Tr::Cell_handle  Cell_handle;

  typedef typename Tds::Edge Etds;
  typedef typename Tds::Edge_iterator  Iterator_base;

  typedef Triangulation_edge_iterator_3<Tr>      Edge_iterator;
public:
  typedef Edge            value_type;
  typedef Edge *          pointer;
  typedef Edge &          reference;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef std::bidirectional_iterator_tag   iterator_category;

  Triangulation_edge_iterator_3()
    : _ib(), _tr(NULL), _inf(true)
  {}
        
  Triangulation_edge_iterator_3(const Tr *tr, bool inf)
    : _ib( &(tr->tds())), _tr(tr), _inf(inf)
  { 
    if (! _inf) {
	while ( // ( _ib != _tr->tds().cells_end() ) &&
	   // useless : there must be at least one finite cell
	   // since precond in _ib : dimension == 3
	       _tr->is_infinite(*_ib))
	  { ++_ib; }
    }
  }
        
  Triangulation_edge_iterator_3(const Tr *tr)
    : _ib( &(tr->tds()), 1), _tr(tr), _inf(true)
    // _inf is initialized but should never be used
  { }
       
  bool
  operator==(const Edge_iterator & ei) const
  {
    if ( _tr != ei._tr ) 
      return false;

    if ( ( _ib == _tr->tds().edges_end() ) 
	 || ( ei._ib == _tr->tds().edges_end() ) ) 
      return ( _ib == ei._ib );

    return ( ( _ib == ei._ib ) && ( _inf == ei._inf ) );
  }

  bool
  operator!=(const Edge_iterator & ei) const
  {
    return !(*this == ei);
  }

  Edge_iterator &
  operator++()
  {
    if (_inf) {
      ++_ib;
    }
    else {
      do {
	++_ib; 
      } while ( _ib != _tr->tds().edges_end() && _tr->is_infinite(*_ib));
    }
    return *this;   
  }

  Edge_iterator &
  operator--()
  {
    if (_inf) {
      --_ib;
    }
    else{
      do {
	--_ib;
      } while ( ( _ib != _tr->tds().edges_end() ) && _tr->is_infinite(*_ib));
    }
    return *this;   
  }

  Edge_iterator
  operator++(int)
  {
    Edge_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Edge_iterator
  operator--(int)
  {
    Edge_iterator tmp(*this);
    --(*this);
    return tmp;
  }
        
  Edge operator*() const
  {
    return *_ib;
  }
     
private:
  Iterator_base _ib;
  const Tr * _tr;
  bool _inf; // if _inf == true, traverses all edges
               // else only traverses finite edges
};

template <class Tr>
class Triangulation_facet_iterator_3
{
  typedef typename Tr::Triangulation_data_structure Tds;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Cell Cell;

  typedef typename Tds::Facet Ftds;
  typedef typename Tds::Facet_iterator  Iterator_base;

  typedef Triangulation_facet_iterator_3<Tr> Facet_iterator;

public:
  typedef Facet           value_type;
  typedef Facet *         pointer;
  typedef Facet &         reference;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef std::bidirectional_iterator_tag   iterator_category;

  Triangulation_facet_iterator_3()
    : _ib(), _tr(NULL), _inf(true)
  {}
        
  Triangulation_facet_iterator_3(const Tr *tr, bool inf)
    : _ib( &(tr->tds())), _tr(tr), _inf(inf)
  {       
    if (! _inf) {
      while ( // ( _ib != _tr->tds().cells_end() ) &&
	     // useless : there must be at least one finite cell
	     // since precond in _ib : dimension == 3
	     _tr->is_infinite(*_ib) )
	{ ++_ib; }
    }
  }

  Triangulation_facet_iterator_3(const Tr *tr)
    : _ib( &(tr->tds()), 1), _tr(tr), _inf(true)
  // _inf is initialized but should never be used
  {}

  bool
  operator==(const Facet_iterator & fi) const
  {
    if ( _tr != fi._tr ) 
      return false;

    if ( ( _ib == _tr->tds().facets_end() ) 
	 || ( fi._ib == _tr->tds().facets_end() ) ) 
      return ( _ib == fi._ib );

    return ( ( _ib == fi._ib ) && ( _inf == fi._inf ) );
  }

  bool
  operator!=(const Facet_iterator & fi) const
  {
    return !(*this == fi);
  }

  Facet_iterator &
  operator++()
  {
    if (_inf) {
      ++_ib;
    }
    else {
      do {
	++_ib; 
      } while ( ( _ib != _tr->tds().facets_end() )
		&& _tr->is_infinite(*_ib) );
    }
    return *this;   
  }

  Facet_iterator &
  operator--()
  {
    if (_inf) {
      --_ib;
    }
    else{
      do {
	--_ib;
      } while ( ( _ib != _tr->tds().facets_end() )
		&& _tr->is_infinite(*_ib) );
    }
    return *this;   
  }

  Facet_iterator
  operator++(int)
  {
    Facet_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Facet_iterator
  operator--(int)
  {
    Facet_iterator tmp(*this);
    --(*this);
    return tmp;
  }
        
  Facet operator*() const
  {
    return *_ib;
  }
     
private:
  Iterator_base   _ib ;
  const Tr * _tr;
  bool _inf; // if _inf == true, traverses all facets
               // else only traverses finite facets
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_ITERATORS_3_H
