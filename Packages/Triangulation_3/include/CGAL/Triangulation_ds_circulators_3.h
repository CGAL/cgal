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
// file          : include/CGAL/Triangulation_ds_circulators_3.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DS_CIRCULATORS_3_H
#define CGAL_TRIANGULATION_DS_CIRCULATORS_3_H

#include <pair.h>
#include <CGAL/triple.h>
#include <CGAL/circulator.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_3.h>


template < class Tds >
class CGAL_Triangulation_ds_cell_circulator_3
  : public CGAL_Bidirectional_circulator_base<typename Tds::Cell, ptrdiff_t, size_t>
{
  // circulates on cells around a given edge
public:

  friend class Tds;
  typedef typename Tds::Cell Cell;
  typedef typename Tds::Edge Edge;
  typedef typename Tds::Vertex Vertex;

  typedef CGAL_Triangulation_ds_cell_circulator_3<Tds> Cell_circulator;

  // CONSTRUCTORS

  CGAL_Triangulation_ds_cell_circulator_3()
    : _Triangulation_data_structure_3(NULL), _e(NULL), pos(NULL), prev(NULL)
  {}
        
  CGAL_Triangulation_ds_cell_circulator_3(Tds * Triangulation_data_structure_3, Edge e)
    : _Triangulation_data_structure_3(Triangulation_data_structure_3), _e(e)
  {
    CGAL_triangulation_precondition( e.first != NULL && 
				     (e.second==0 || e.second==1 ||
				      e.second==2 || e.second==3 ) &&
				     (e.third==0 || e.third==1 ||
				      e.third==2 || e.third==3 ) );
//     if ( _Triangulation_data_structure_3->dimension() <3 ) useless since precondition in Triangulation_data_structure_3 incident_cells
//       {
// 	pos = NULL;
// 	prev = NULL;
//       }
//     else {
      pos = e.first;
      int k = other(e.second, e.third);
      prev = e.first->neighbor(k);
//    }
  }

  CGAL_Triangulation_ds_cell_circulator_3(const Cell_circulator & ccir)
    : _Triangulation_data_structure_3(ccir._Triangulation_data_structure_3), _e(ccir._e), pos(ccir.pos), prev(ccir.prev)
  {}
        
  Cell_circulator & operator++()
  {
    CGAL_triangulation_precondition( (pos != NULL) );
    //then dimension() cannot be < 3
    int i= pos->index(_e.first->vertex(_e.second));
    int j = pos->index(_e.first->vertex(_e.third));
    int k = other(i,j);
    Cell * tmp = pos;
    Cell * n = pos->neighbor(k);
    if ( n == prev ) {
      pos = pos->neighbor( 6-i-j-k );
    }
    else pos = n;

    prev = tmp;
    return *this;
  }
        
  Cell_circulator operator++(int)
  {
    CGAL_triangulation_precondition( (pos  != NULL) );
    Cell_circulator tmp(*this);
    ++(*this);
    return tmp;
  }
        
  Cell_circulator& operator--()
  {
    CGAL_triangulation_precondition( (pos != NULL) ); // then prev != NULL too
    int i = prev->index(_e.first->vertex(e.second));
    int j = prev->index(_e.first->vertex(e.third));
    int k = other(i,j);
    Cell * tmp = prev;
    Cell * n = prev->neighbor(k);
    if ( n == pos ) {
      prev = prev->neighbor( 6-i-j-k );
    }
    else prev = n;
    pos = tmp;
    return *this;
  }
        
        
  Cell_circulator operator--(int)
  {
    CGAL_triangulation_precondition( (pos != NULL) && (_v != NULL) );
    Cell_circulator tmp(*this);
    --(*this);
    return tmp;
  }

  inline Cell& operator*() const
  {
    return *pos;
  }

  inline Cell* operator->() const
  {
    return pos;
  }
  
  bool operator==(const Cell_circulator & ccir) const
  {
    return ( (_Triangulation_data_structure_3 == ccir._Triangulation_data_structure_3) && (_e == ccir._e) && 
	     (pos == ccir.pos) && (prev == ccir.prev) );
  }
        
        
  bool operator!=(const Cell_circulator & ccir) const
  {
    return ! (*this == ccir);
  }
        
//   bool
//   operator==(CGAL_NULL_TYPE n) const
//   {
//     CGAL_triangulation_assertion( n == NULL);
//     return (pos == NULL);
//   }
        
//   bool
//   operator!=(CGAL_NULL_TYPE n) const
//   {
//     CGAL_triangulation_assertion( n == NULL);
//     return ! (*this == NULL);
//   }
        
private: 
  Tds* _Triangulation_data_structure_3;
  Edge _e;
  Cell* pos; // current cell
  Cell* prev; // preceding cell

  // given i and j supposed to be different and in 0,1,2,3
  // i supposed to be < j
  // computes the index of a vertex different from i and j in a cell 
  // works in dimension 2 or 3
  static int other(int i, int j)
  {
    int min, max;
    if (i<j) {
      min = i;
      max = j;
    }
    else {
      min = j;
      max = i;
    }
    if ( min == 0 ) {
      if ( max == 1 ) return 2;
      else return 1;
    }
    else return 0;
  }
};

#endif CGAL_TRIANGULATION_DS_CIRCULATORS_3_H

