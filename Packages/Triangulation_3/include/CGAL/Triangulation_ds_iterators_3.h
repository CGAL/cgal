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
// file          : include/CGAL/Triangulation_ds_iterators_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DS_ITERATORS_3_H
#define CGAL_TRIANGULATION_DS_ITERATORS_3_H

#include <pair.h>
#include <CGAL/triple.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_3.h>

#include <CGAL/Triangulation_ds_circulators_3.h>

template < class Vb, class Cb >
class CGAL_Triangulation_ds_vertex_3;
template < class Vb, class Cb >
class CGAL_Triangulation_ds_cell_3;
template < class Vb, class Cb >
class CGAL_Triangulation_data_structure_3;
template <class Tds>
class CGAL_Triangulation_ds_cell_circulator_3;

template<class Tds>
class CGAL_Triangulation_ds_cell_iterator_3
  : public bidirectional_iterator<typename Tds::Cell, ptrdiff_t>
{
public:
  typedef typename Tds::Cell  Cell;

  typedef CGAL_Triangulation_ds_cell_iterator_3<Tds> Cell_iterator;

  // CONSTRUCTORS

  CGAL_Triangulation_ds_cell_iterator_3()
    : _Triangulation_data_structure_3(NULL), pos(NULL)
    {}
  
  CGAL_Triangulation_ds_cell_iterator_3(Tds * Triangulation_data_structure_3)
    : _Triangulation_data_structure_3(Triangulation_data_structure_3)
    {
      if ( _Triangulation_data_structure_3->dimension() <3 ) { pos = NULL; } // there is no cell yet
      else { pos = _Triangulation_data_structure_3->list_of_cells()._next_cell; } // first cell of the list
    }

  // used to initialize the past-the end iterator
  CGAL_Triangulation_ds_cell_iterator_3(Tds* Triangulation_data_structure_3, int i)
    : _Triangulation_data_structure_3(Triangulation_data_structure_3)
    {
      if ( _Triangulation_data_structure_3->dimension() <3 ) { pos = NULL; } // there is no cell yet
      else { pos = _Triangulation_data_structure_3->past_end_cell(); }
    }
  
// MT ne sert a rien ?
//   CGAL_Triangulation_ds_iterator_base_3(Tds* Triangulation_data_structure_3, Cell* f)
//     : pos(f), _Triangulation_data_structure_3(Triangulation_data_structure_3)
//     {}
//   CGAL_Triangulation_ds_cell_iterator_3(const Cell_iterator& fi)
//     : Iterator_base(fi._Triangulation_data_structure_3, fi.pos)
//     {}
  
  // OPERATORS
  
  Cell_iterator &
  operator=(const Cell_iterator & fi)
    {
      pos = fi.pos;
      _Triangulation_data_structure_3 = fi._Triangulation_data_structure_3;
      return *this;
    }

  bool
  operator==(const Cell_iterator & fi) const
    {
      return ((pos == fi.pos )&&(_Triangulation_data_structure_3==fi._Triangulation_data_structure_3));
    }
  
  bool
  operator!=(const Cell_iterator& fi)
    {
      return !(*this == fi);
    }
  
  Cell_iterator &
  operator++()
    {
      // the user should not use ++ on a non-initialized iterator
//       if (pos==NULL){
// 	return *this;            
//       }
      // removed : past the end iterator can go further
//       if ( pos == _Triangulation_data_structure_3->past_end_cell() ){
// 	return *this;    
//       } 
      pos = pos->_next_cell;
      return *this;
    }

  Cell_iterator &
  operator--()
    {
      // the user should not use ++ on a non-initialized iterator
//       if (pos==NULL){
// 	return *this;            
//       }
      // removed : begin iterator can go backwards
//       if ( pos == _Triangulation_data_structure_3->list_of_cells()._next_cell ) {
// 	return *this;
//       }
      pos = pos->_previous_cell;
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
        
  inline Cell & operator*() const
    {
      return *pos;
    }
    
  inline Cell* operator->() const
    {
      return pos;
    }

private:
  Tds*  _Triangulation_data_structure_3;
  Cell* pos;

};

template < class Tds>
class CGAL_Triangulation_ds_vertex_iterator_3
  : public bidirectional_iterator<typename Tds::Vertex, ptrdiff_t>
{
// traverses the list of cells and report for each cell 
// the vertices whose cell() is the current cell

public:

  typedef typename Tds::Vertex Vertex;
  typedef typename Tds::Cell  Cell;
  typedef CGAL_Triangulation_ds_vertex_iterator_3<Tds> Vertex_iterator;
  
  CGAL_Triangulation_ds_vertex_iterator_3()
    : _Triangulation_data_structure_3(NULL), pos(NULL), index(0)
    {}
  
  CGAL_Triangulation_ds_vertex_iterator_3(Tds * Triangulation_data_structure_3)
    : _Triangulation_data_structure_3(Triangulation_data_structure_3), index(0)
    {
      if ( _Triangulation_data_structure_3->number_of_vertices() == 0 ) { pos = NULL; }
      else { 
	pos = Triangulation_data_structure_3->list_of_cells()._next_cell; 
	while ( (pos != _Triangulation_data_structure_3->past_end_cell())
		  && (pos->vertex(index)->cell() != pos) ) {
	  increment();
	}
      }
    }
  
  // used to initialize the past-the end iterator
  CGAL_Triangulation_ds_vertex_iterator_3(Tds* Triangulation_data_structure_3, int i)
    : _Triangulation_data_structure_3(Triangulation_data_structure_3), index(0)
    {
      if ( _Triangulation_data_structure_3->number_of_vertices() == 0 ) { pos = NULL; }
      else { pos = Triangulation_data_structure_3->past_end_cell() ; }
    }
  
  Vertex_iterator&
  operator++()
  {
    //    if (_Triangulation_data_structure_3->number_of_vertices() == 0){ return *this;} should not be accessed
    do {
      increment();
    } while ( (pos != _Triangulation_data_structure_3->past_end_cell())
	      && (pos->vertex(index)->cell() != pos) );
    return *this;
  }
    
  Vertex_iterator&
  operator--()
  {
    // if (_Triangulation_data_structure_3->number_of_vertices() == 0) { return *this;} should not be accessed
    do{
      if (index == 0){
	// all the vertices of the current cell have been examined
	int d = _Triangulation_data_structure_3->dimension();
	if ( d >= 0 ) {index = d;}
	else {index = 0;}
	pos = pos->_previous_cell;
      }
      else { index--; }
    } while ( (pos != _Triangulation_data_structure_3->past_end_cell())
	      && (pos->vertex(index)->cell() != pos) );
    return *this;
  }
    
  Vertex_iterator operator++(int)
    {
      Vertex_iterator tmp(*this);
      ++(*this);
      return tmp;
    }
    
  Vertex_iterator operator--(int)
    {
      Vertex_iterator tmp(*this);
      --(*this);
      return tmp;
    }
    
  bool operator==(const Vertex_iterator& fi) const
    {
      return ( (_Triangulation_data_structure_3 == fi._Triangulation_data_structure_3) && (pos == fi.pos) && (index == fi.index) );
    }
    
  bool operator!=(const Vertex_iterator& fi) const
    {
      return !(*this == fi);
    }
    
  inline Vertex & operator*() const
    {
      // case pos == NULL should not be accessed, there is no vertex
      return *(pos->vertex(index));
    }
    
  inline Vertex* operator->() const
    {
      // case pos == NULL should not be accessed, there is no vertex
      return pos->vertex(index);
    }

private:
  Tds*  _Triangulation_data_structure_3;
  Cell* pos; // current "cell". Even if the dimension is <3 when 
              // there is no true cell yet.
  int index; // index of the current vertex in the current cell

  void 
  increment()
  {
    if (index >= _Triangulation_data_structure_3->dimension()) {
      // all the vertices of the current cell have been examined
      index = 0;
      pos = pos->_next_cell;
    }
    // be careful : index should always be 0 when pos = past_end_cell
    else { index++; }
  }

};

template < class Tds>
class CGAL_Triangulation_ds_facet_iterator_3
  : public bidirectional_iterator<typename Tds::Facet, ptrdiff_t>
{
// traverses the list of cells and report for each cell 
// the vertices whose cell() is the current cell

public:

  typedef typename Tds::Cell Cell;
  typedef typename Tds::Facet Facet;
  typedef CGAL_Triangulation_ds_facet_iterator_3<Tds> Facet_iterator;
  
  CGAL_Triangulation_ds_facet_iterator_3()
    : _Triangulation_data_structure_3(NULL), pos(NULL), index(0)
    {}
  
  CGAL_Triangulation_ds_facet_iterator_3(Tds * Triangulation_data_structure_3)
    : _Triangulation_data_structure_3(Triangulation_data_structure_3), index(0)
    {
      switch ( _Triangulation_data_structure_3->dimension() ) {
      case 2:
	pos = Triangulation_data_structure_3->list_of_cells()._next_cell; 
	index = 3;
	return;
      case 3:
	pos = Triangulation_data_structure_3->list_of_cells()._next_cell; 
	while ( // useless (pos != _Triangulation_data_structure_3->past_end_cell()) &&
	       // there must be at least one facet
	       ( pos->neighbor(index) < pos ) ) {
	  increment();
	}
	return;
      default:
	pos = NULL;
	return;
      }
    }
  
  // used to initialize the past-the end iterator
  CGAL_Triangulation_ds_facet_iterator_3(Tds* Triangulation_data_structure_3, int i)
    : _Triangulation_data_structure_3(Triangulation_data_structure_3), index(0)
    {
      if ( _Triangulation_data_structure_3->dimension() < 2 ) { pos = NULL; }
      else { 
	pos = Triangulation_data_structure_3->past_end_cell() ; 
	if  ( _Triangulation_data_structure_3->dimension() == 2 ) { index = 3; }
      }
    }
  
  Facet_iterator&
  operator++()
  {
    //    if (_Triangulation_data_structure_3->dimension() < 2){ return *this;} should not be accessed
    if ( _Triangulation_data_structure_3->dimension() == 3 ) {
      do {
	increment();
      } while ( (pos != _Triangulation_data_structure_3->past_end_cell())
		&& ( (pos->neighbor(index)) < pos ) );
      // reports a facet when the current cell has a pointer inferior
      // to the pointer of the neighbor cell
    }
    else { pos = pos->_next_cell; }
    return *this;
  }
    
  Facet_iterator&
  operator--()
  {
    // if (_Triangulation_data_structure_3->dimension() < 2) { return *this;} should not be accessed
    if ( _Triangulation_data_structure_3->dimension() == 2 ) {
      pos = pos->_previous_cell; // index remains 3
    }
    else { // dimension 3
      do{
	if (index == 0){
	  // all the facets of the current cell have been examined
	  index = 3;
	  pos = pos->_previous_cell;
	}
	else { index--; }
      } while ( (pos != _Triangulation_data_structure_3->past_end_cell())
		&& (pos->neighbor(index) < pos ) );
      // reports a facet when the current cell has a pointer inferior
      // to the pointer of the neighbor cell
    }
    return *this;
  }
    
  Facet_iterator operator++(int)
    {
      Facet_iterator tmp(*this);
      ++(*this);
      return tmp;
    }
    
  Facet_iterator operator--(int)
    {
      Facet_iterator tmp(*this);
      --(*this);
      return tmp;
    }
    
  bool operator==(const Facet_iterator& fi) const
    {
      return ( (_Triangulation_data_structure_3 == fi._Triangulation_data_structure_3) && (pos == fi.pos) && (index == fi.index) );
    }
    
  bool operator!=(const Facet_iterator& fi) const
    {
      return !(*this == fi);
    }
    
  Facet
  operator*() const
    {
      // case pos == NULL should not be accessed, there is no facet when dimension <2
      return make_pair(pos, index);
    }
    
private:
  Tds*  _Triangulation_data_structure_3;
  Cell* pos; // current "cell". Even if the dimension is <3 when 
              // there is no true cell yet.
  int index; // index of the current facet in the current cell

  void
  increment()
  {
//     if ( _Triangulation_data_structure_3->dimension() == 2 ) {
//       pos = pos->_next_cell; // index remains 3
//     }
//     else { // dimension 3
    if (index == 3) {
      // all the feces of the current cell have been examined
      index = 0;
      pos = pos->_next_cell;
    }
    // be careful : index should always be 0 when pos = past_end_cell
    else { index++; }
//    }
  }

};

template < class Tds>
class CGAL_Triangulation_ds_edge_iterator_3
  : public bidirectional_iterator<typename Tds::Edge, ptrdiff_t>
{
// traverses the list of cells and report for each cell 
// the vertices whose cell() is the current cell

public:

  typedef typename Tds::Cell Cell;
  typedef typename Tds::Edge Edge;
  typedef CGAL_Triangulation_ds_edge_iterator_3<Tds> Edge_iterator;
  typedef CGAL_Triangulation_ds_cell_circulator_3<Tds> Cell_circulator;
  
  CGAL_Triangulation_ds_edge_iterator_3()
    : _Triangulation_data_structure_3(NULL), pos(NULL), b(0), e(1)
    {}
  
  CGAL_Triangulation_ds_edge_iterator_3(Tds * Triangulation_data_structure_3)
    : _Triangulation_data_structure_3(Triangulation_data_structure_3), b(0), e(1)
    {
      switch ( _Triangulation_data_structure_3->dimension() ) {
      case 1:
	{
	  pos = Triangulation_data_structure_3->list_of_cells()._next_cell;
	  return;
	}
      case 2:
	{
	  pos = Triangulation_data_structure_3->list_of_cells()._next_cell; 
	  while ( // useless (pos != _Triangulation_data_structure_3->past_end_cell()) && 
		 // there must be at least one edge
		 ( pos->neighbor(3-b-e) < pos) ) {
	    increment2();
	  }
	  return;
	}
      case 3:
	{
	  pos = Triangulation_data_structure_3->list_of_cells()._next_cell;
	  bool notfound = true;
	  while ( // useless (pos != _Triangulation_data_structure_3->past_end_cell()) &&
		 // there must be at least one edge
		 notfound ) {
	    Cell_circulator ccir = _Triangulation_data_structure_3->incident_cells(CGAL_make_triple(pos,b,e));
	    do {
	      ++ccir;
	    } while ( &(*ccir) > pos ); 
	    // loop terminates since it stops at least when ccir = pos
	    if ( &(*ccir) == pos ) {// pos is the cell with minimum pointer
	      notfound = false;
	    }
	    else {
	      increment3();
	    }
	  }
	  return;
	}
      default:
	{
	  pos = NULL;
	  return;
	}
      }
    }
  
  // used to initialize the past-the end iterator
  CGAL_Triangulation_ds_edge_iterator_3(Tds* Triangulation_data_structure_3, int i)
    : _Triangulation_data_structure_3(Triangulation_data_structure_3), b(0), e(1)
    {
      if ( _Triangulation_data_structure_3->dimension() < 1 ) { pos = NULL; }
      else { 
	pos = Triangulation_data_structure_3->past_end_cell() ; 
	// always b == 0 and e == 1 for past-end-cell
// 	switch ( _Triangulation_data_structure_3->dimension() ) {
// 	  // b and e set to the indices of the "last" --the
// 	  // edges are lexicographically ordered-- possible edge in a cell
// 	case 1:
// 	  b = 0; 
// 	  e = 1;
// 	  break;
// 	case 2:
// 	  b = 0; //2
// 	  e = 1; //0
// 	  break;
// 	case 3:
// 	  b = 0; //2
// 	  e = 1; //3
// 	  break;
// 	}
      }
    }
  
  Edge_iterator&
  operator++()
  {
    switch ( _Triangulation_data_structure_3->dimension() ) {
      // dimension < 1 should not occur
    case 1:
      pos = pos->_next_cell;
      break;
    case 2:
      do {
	increment2();
      } while ( (pos != _Triangulation_data_structure_3->past_end_cell()) && 
		( pos->neighbor(3-b-e) < pos) );
      break;
    case 3:
      bool notfound = true;
      do {
	increment3();
	if (pos != _Triangulation_data_structure_3->past_end_cell()) {
	  Cell_circulator ccir = _Triangulation_data_structure_3->incident_cells(CGAL_make_triple(pos,b,e));
	  do {
	    ++ccir;
	  } while ( &(*ccir) > pos );
	  if ( &(*ccir) == pos ) {// pos is the cell with minimum pointer
	    notfound = false;
	  }
	}
      } while ( (pos != _Triangulation_data_structure_3->past_end_cell()) &&
		notfound );
      break;
    }
    return *this;
  }
    
  Edge_iterator&
  operator--()
  {
    // if (_Triangulation_data_structure_3->dimension() < 2) { return *this;} should not be accessed
    switch (  _Triangulation_data_structure_3->dimension() ) {
    case 1:
      pos = pos->_previous_cell; // b, e remain 0, 1
      break;
    case 2:
      do {
	if (b == 0) {
	  b = 2; e = 0;
	  pos = pos->_previous_cell;
	}
	else {
	  b--; 
	  e = b+1; // case b==2, e==0 forbids to write e--
	}
      } while ( (pos != _Triangulation_data_structure_3->past_end_cell()) && 
		( pos->neighbor(3-b-e) < pos) );
      break;
    case 3:
      bool notfound = true;
      do {
	if (b == 0) {
	  if (e == 1) {
	    // all the edges of the current cell have been examined
	    b = 2; e = 3;
	    pos = pos->_previous_cell;
	  }
	  else {
	    e--;
	  }
	}
	else {
	  if (e == b+1) {
	    b--;
	    e = 3;
	  }
	  else {
	    e--;
	  }
	}
	Cell_circulator ccir = _Triangulation_data_structure_3->incident_cells(CGAL_make_triple(pos,b,e));
	while ( &(*ccir) > pos ) {
	  --ccir;
	}
	if ( &(*ccir) == pos ) {// pos is the cell with minimum pointer
	  notfound = false;
	}
      } while ( (pos != _Triangulation_data_structure_3->past_end_cell()) &&
		notfound );
      break;
    }
    // reports an edge when the current cell has a pointer inferior
    // to the pointer of the neighbor cell
    return *this;
  }
    
  Edge_iterator operator++(int)
    {
      Edge_iterator tmp(*this);
      ++(*this);
      return tmp;
    }
    
  Edge_iterator operator--(int)
    {
      Edge_iterator tmp(*this);
      --(*this);
      return tmp;
    }
    
  bool operator==(const Edge_iterator& ei) const
    {
      return ( (_Triangulation_data_structure_3 == ei._Triangulation_data_structure_3) && (pos == ei.pos) 
	       && (b == ei.b) && (e == ei.e) );
    }
    
  bool operator!=(const Edge_iterator& ei) const
    {
      return !(*this == ei);
    }
    
  Edge
  operator*() const
    {
      // case pos == NULL should not be accessed, there is no edge when dimension <2
      return CGAL_make_triple(pos, b, e);
    }
    
private:
  Tds*  _Triangulation_data_structure_3;
  Cell* pos; // current "cell". Even if the dimension is <3 when 
              // there is no true cell yet.
  int b; // index of the first endpoint of the current edge in the current cell
  int e; // index of the second endpoint of the current edge in the current cell

  void
  increment2()
  {
//     switch ( _Triangulation_data_structure_3->dimension() ) {
//     case 1:
//       pos = pos->_next_cell; // b, e remain 0, 1
//     case 2: { 
    if (b == 2) { // e == 0
      // all the edges of the current cell have been examined
      b = 0; e = 1;
      pos = pos->_next_cell;
    }
    // be careful : index should always be 0 when pos = past_end_cell
    else { 
      b++; 
      if ( b == 2 ) {
	e = 0;
      }
      else { // b==1
	e = 2;
      }
    }
//    }
  }

  void
  increment3()
  {
//    case 3: { 
    if (b == 2) { // then e == 3
      // all the edges of the current cell have been examined
      b = 0; e = 1;
      pos = pos->_next_cell;
    }
    else {
      if (e == 3) {
	b++;
	e = b+1;
      }
      else {
	e++;
      }
    }
//    }
//    }
  }

};
#endif CGAL_TRIANGULATION_DS_ITERATORS_3_H

