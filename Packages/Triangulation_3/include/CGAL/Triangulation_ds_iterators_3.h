// ============================================================================
//
// $Id$
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DS_ITERATORS_H
#define CGAL_TRIANGULATION_DS_ITERATORS_H

#include <pair.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names.h>

template < class Vb, class Cb >
class CGAL_Triangulation_ds_vertex;

template < class Vb, class Cb >
class CGAL_Triangulation_ds_cell;

template < class Vb, class Cb >
class CGAL_Triangulation_data_structure;


template<class Tds>
class CGAL_Triangulation_ds_cell_iterator
  : public bidirectional_iterator<typename Tds::Cell, ptrdiff_t>
{
public:
  typedef typename Tds::Cell  Cell;

  typedef CGAL_Triangulation_ds_cell_iterator<Tds> Cell_iterator;

  // CONSTRUCTORS

  CGAL_Triangulation_ds_cell_iterator()
    : pos(NULL), _tds(NULL)
    {}
  
  CGAL_Triangulation_ds_cell_iterator(Tds * tds)
    : _tds(tds)
    {
      if ( _tds->dimension() <3 ) { pos = NULL; } // there is no cell yet
      else { pos = _tds->list_of_cells()._next_cell; } // first cell of the list
    }

  // used to initialize the past-the end iterator
  CGAL_Triangulation_ds_cell_iterator(Tds* tds, int i)
    : _tds(tds)
    {
      if ( _tds->dimension() <3 ) { pos = NULL; } // there is no cell yet
      else { pos = _tds->past_end_cell(); }
    }
  
// MT ne sert a rien ?
//   CGAL_Triangulation_ds_iterator_base(Tds* tds, Cell* f)
//     : pos(f), _tds(tds)
//     {}
//   CGAL_Triangulation_ds_cell_iterator(const Cell_iterator& fi)
//     : Iterator_base(fi._tds, fi.pos)
//     {}
  
  // OPERATORS
  
  Cell_iterator &
  operator=(const Cell_iterator & fi)
    {
      pos = fi.pos;
      _tds = fi._tds;
      return *this;
    }

  bool
  operator==(const Cell_iterator & fi) const
    {
      return ((pos == fi.pos )&&(_tds==fi._tds));
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
//       if ( pos == _tds->past_end_cell() ){
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
//       if ( pos == _tds->list_of_cells()._next_cell ) {
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
  Tds*  _tds;
  Cell* pos;

};

template < class Tds>
class CGAL_Triangulation_ds_vertex_iterator
  : public bidirectional_iterator<typename Tds::Vertex, ptrdiff_t>
{
// traverses the list of cells and report for each cell 
// the vertices whose cell() is the current cell

public:

  typedef typename Tds::Vertex Vertex;
  typedef typename Tds::Cell  Cell;
  typedef CGAL_Triangulation_ds_vertex_iterator<Tds> Vertex_iterator;
  
  CGAL_Triangulation_ds_vertex_iterator()
    : _tds(NULL), pos(NULL), index(0)
    {}
  
  CGAL_Triangulation_ds_vertex_iterator(Tds * tds)
    : _tds(tds), index(0)
    {
      if ( _tds->number_of_vertices() == 0 ) { pos = NULL; }
      else { 
	pos = tds->list_of_cells()._next_cell; 
	while ( (pos != _tds->past_end_cell())
		  && (pos->vertex(index)->cell() != pos) ) {
	  increment();
	}
      }
    }
  
  // used to initialize the past-the end iterator
  CGAL_Triangulation_ds_vertex_iterator(Tds* tds, int i)
    : _tds(tds), index(0)
    {
      if ( _tds->number_of_vertices() == 0 ) { pos = NULL; }
      else { pos = tds->past_end_cell() ; }
    }
  
  Vertex_iterator&
  operator++()
  {
    //    if (_tds->number_of_vertices() == 0){ return *this;} should not be accessed
    do {
      increment();
    } while ( (pos != _tds->past_end_cell())
	      && (pos->vertex(index)->cell() != pos) );
    return *this;
  }
    
  Vertex_iterator&
  operator--()
  {
    // if (_tds->number_of_vertices() == 0) { return *this;} should not be accessed
    do{
      if (index == 0){
	// all the vertices of the current cell have been examined
	int d = _tds->dimension();
	if ( d >= 0 ) {index = d;}
	else {index = 0;}
	pos = pos->_previous_cell;
      }
      else { index--; }
    } while ( (pos != _tds->past_end_cell())
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
      return ( (_tds == fi._tds) && (pos == fi.pos) && (index == fi.index) );
    }
    
  bool operator!=(const Vertex_iterator& fi) const
    {
      return !(*this == fi);
    }
    
  inline Vertex& operator*() const
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
  Tds*  _tds;
  Cell* pos; // current "cell". Even if the dimension is <3 when 
              // there is no true cell yet.
  int index; // index of the current vertex in the current cell

  void 
  increment()
  {
    if (index >= _tds->dimension()) {
      // all the vertices of the current cell have been examined
      index = 0;
      pos = pos->_next_cell;
    }
    // be careful : index should always be 0 when pos = past_end_cell
    else { index++; }
  }

};

template < class Tds>
class CGAL_Triangulation_ds_facet_iterator
  : public bidirectional_iterator<typename Tds::Facet, ptrdiff_t>
{
// traverses the list of cells and report for each cell 
// the vertices whose cell() is the current cell

public:

  typedef typename Tds::Cell Cell;
  typedef typename Tds::Facet Facet;
  typedef CGAL_Triangulation_ds_facet_iterator<Tds> Facet_iterator;
  
  CGAL_Triangulation_ds_facet_iterator()
    : _tds(NULL), pos(NULL), index(0)
    {}
  
  CGAL_Triangulation_ds_facet_iterator(Tds * tds)
    : _tds(tds), index(0)
    {
      switch ( _tds->dimension() ) {
      case 2:
	pos = tds->list_of_cells()._next_cell; 
	index = 3;
	return;
      case 3:
	while ( (pos != _tds->past_end_cell())
		&& ( (pos->neighbor(index)) < pos ) ) {
	  increment();
	}
	return;
      default:
	pos = NULL;
	return;
      }
    }
  
  // used to initialize the past-the end iterator
  CGAL_Triangulation_ds_facet_iterator(Tds* tds, int i)
    : _tds(tds), index(0)
    {
      if ( _tds->dimension() < 2 ) { pos = NULL; }
      else { 
	pos = tds->past_end_cell() ; 
	if  ( _tds->dimension() == 2 ) { index = 3; }
      }
    }
  
  Facet_iterator&
  operator++()
  {
    //    if (_tds->dimension() < 2){ return *this;} should not be accessed
    if ( _tds->dimension() == 3 ) {
      do {
	increment();
      } while ( (pos != _tds->past_end_cell())
		&& ( (pos->neighbor(index)) < pos ) );
      // reports a facet when the current cell has a pointer superior 
      // to the pointer of the neighbor cell
    }
    else { increment(); }
    return *this;
  }
    
  Facet_iterator&
  operator--()
  {
    // if (_tds->dimension() < 2) { return *this;} should not be accessed
    do{
      if ( _tds->dimension() == 2 ) {
	pos = pos->_previous_cell; // index remains 3
      }
      else { // dimension 3
	if (index == 0){
	  // all the facets of the current cell have been examined
	  index = 3;
	  pos = pos->_previous_cell;
	}
	else { index--; }
      }
    } while ( (pos != _tds->past_end_cell())
	      && (pos->neighbor(index) < pos ) );
    // reports a facet when the current cell has a pointer superior 
    // to the pointer of the neighbor cell
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
      return ( (_tds == fi._tds) && (pos == fi.pos) && (index == fi.index) );
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
  Tds*  _tds;
  Cell* pos; // current "cell". Even if the dimension is <3 when 
              // there is no true cell yet.
  int index; // index of the current facet in the current cell

  void
  increment()
  {
    if ( _tds->dimension() == 2 ) {
      pos = pos->_next_cell; // index remains 3
    }
    else { // dimension 3
      if (index == 3) {
	// all the feces of the current cell have been examined
	index = 0;
	pos = pos->_next_cell;
      }
      // be careful : index should always be 0 when pos = past_end_cell
      else { index++; }
    }
  }

};

#endif CGAL_TRIANGULATION_DS_ITERATORS_H

