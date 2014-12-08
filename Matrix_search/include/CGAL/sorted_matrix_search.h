// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#ifndef CGAL_SORTED_MATRIX_SEARCH_H
#define CGAL_SORTED_MATRIX_SEARCH_H 1

#include <CGAL/basic.h>
#include <CGAL/Optimisation/assertions.h>
#include <boost/bind.hpp>
#include <algorithm>
#include <vector>
#include <CGAL/Sorted_matrix_search_traits_adaptor.h>

namespace CGAL {
template < class Matrix >
class Padded_matrix {
public:
  typedef typename Matrix::Value Value;

  Padded_matrix( const Matrix& m) : matrix( &m) {}

  Value
  operator()( int x, int y) const
  // padded access operator
  {
    return matrix->operator()(
      x < matrix->number_of_columns() ?
        x : matrix->number_of_columns() - 1,
      y < matrix->number_of_rows() ?
        y : matrix->number_of_rows() - 1);
  }

  bool
  is_sorted()
  // tests iff this matrix is sorted, i.e. in each column/row
  // the elements appear in increasing order
  // time complexity is proportional to the number of elements
  {
    for ( int i = 0; i < matrix->number_of_columns(); ++i)
      for ( int j = 0; j < matrix->number_of_rows(); ++j) {
        if ( i > 0 && (*matrix)( i - 1, j) > (*matrix)( i, j))
          return false;
        if ( j > 0 && (*matrix)( i, j - 1) > (*matrix)( i, j))
          return false;
      }
    return true;
  }

private:
  const Matrix* matrix;
};
template < class PaddedMatrix >
class Matrix_cell {
public:
  typedef typename PaddedMatrix::Value Value;

  Matrix_cell(PaddedMatrix m, int xpos = 0, int ypos = 0)
  : base_matrix(m), x(xpos), y(ypos)
  {}

  Value
  min
  BOOST_PREVENT_MACRO_SUBSTITUTION
  () const
  { return base_matrix(x, y); }

  Value
  max
  BOOST_PREVENT_MACRO_SUBSTITUTION
  (int offset) const
  // offset denotes the cell's dimension
  { return base_matrix(x + offset - 1, y + offset - 1); }

  int          x_min() const  { return x; }
  int          y_min() const  { return y; }
  PaddedMatrix matrix() const { return base_matrix; }

  void
  output(std::ostream& o, int dim) const
  {
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < dim; ++j)
        o << base_matrix(x + i, y + j) << " ";
      o << std::endl;
    }
  }

  bool
  check_for(Value v, int dim) const {
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) {
        if (CGAL_NTS abs(base_matrix(x + i, y + j) - v) < Value(1E-10))
          std::cerr << "***" << base_matrix(x + i, y + j) << std::endl;
        if (base_matrix(x + i, y + j) == v)
          return true;
      }
    return false;
  }

private:
  PaddedMatrix base_matrix;
  int x;
  int y;
};
template < class Cell >
struct Cell_min
: public std::unary_function< Cell, typename Cell::Value >
{
  typename Cell::Value
  operator()( const Cell& c) const
  { return (c.min)(); }
};

template < class Cell >
struct Cell_max
: public std::unary_function< Cell, typename Cell::Value > {

  Cell_max( int offset) : ofs( offset) {}

  typename Cell::Value
  operator()( const Cell& c) const
  { return (c.max)( ofs); }

private:
  int ofs;
};


template < class InputIterator, class Traits >
typename Traits::Value
sorted_matrix_search(InputIterator f, InputIterator l, Traits t)
{
  BOOST_USING_STD_MAX();
  using std::iter_swap;
  using std::find_if;
  using std::remove_if;
  using std::logical_or;
  using std::equal_to;
  
  typedef typename Traits::Matrix                   Matrix;
  typedef typename Traits::Value                    Value;
  typedef Padded_matrix< Matrix >                   PaddedMatrix;
  typedef Matrix_cell< PaddedMatrix >               Cell;
  typedef std::vector< Cell >                       Cell_container;
  typedef typename Cell_container::iterator         Cell_iterator;
  // typedef typename Cell_container::reverse_iterator Cell_reverse_iterator;
  
  Cell_container active_cells;
  
  // set of input matrices must not be empty:
  CGAL_optimisation_precondition( f != l);
  
  // for each input matrix insert a cell into active_cells:
  InputIterator i( f);
  int maxdim( -1);
  while ( i != l) {
    CGAL_optimisation_expensive_precondition(
      PaddedMatrix( *i).is_sorted());
    active_cells.push_back( Cell( PaddedMatrix( *i)));
    maxdim = max BOOST_PREVENT_MACRO_SUBSTITUTION ( max BOOST_PREVENT_MACRO_SUBSTITUTION ( (*i).number_of_columns(),
                       (*i).number_of_rows()),
                  maxdim);
    ++i;
  }
  CGAL_optimisation_precondition( maxdim > 0);
  
  
  // current cell dimension:
  int ccd( 1);
  // set ccd to a power of two >= maxdim:
  while ( ccd < maxdim)
    ccd <<= 1;
  
  
  

  // now start the search:

  for (;;) {
  
    if ( ccd > 1) {
      // ------------------------------------------------------
      // divide cells:
      ccd >>= 1;
    
    
      // reserve is required here!
      // otherwise one of the insert operations might cause
      // a reallocation invalidating j
      // (should typically result in a segfault)
      // active_cells.reserve( 4 * active_cells.size());
  
      for ( int j = static_cast<int>(active_cells.size()) - 1 ; j >= 0 ; -- j )
      {
      //for ( Cell_reverse_iterator j( active_cells.rbegin());
      //      j != active_cells.rend();
      //      ++j) {

        Cell lRefCell = active_cells.at(j) ;
        
        // upper-left quarter:
        // Cell( (*j).matrix(),
        //       (*j).x_min(),
        //       (*j).y_min()) remains in active_cells,
        // since it is implicitly shortened by decreasing ccd
    
        // lower-left quarter:
        active_cells.push_back(
          Cell( lRefCell.matrix(),
                lRefCell.x_min(),
                lRefCell.y_min() + ccd));
    
        // upper-right quarter:
        active_cells.push_back(
          Cell( lRefCell.matrix(),
                lRefCell.x_min() + ccd,
                lRefCell.y_min()));
    
        // lower-right quarter:
        active_cells.push_back(
          Cell( lRefCell.matrix(),
                lRefCell.x_min() + ccd,
                lRefCell.y_min() + ccd));
    
      } // for all active cells
    } // if ( ccd > 1)
    else if ( active_cells.size() <= 1) //!!! maybe handle <= 3
      break;
      
    // there has to be at least one cell left:
    CGAL_optimisation_assertion( active_cells.size() > 0);
    
    // ------------------------------------------------------
    // compute medians of smallest and largest elements:
    
    
    int lower_median_rank = static_cast<int>((active_cells.size() - 1) >> 1);
    int upper_median_rank = static_cast<int>(active_cells.size() >> 1);
    

    // compute upper median of cell's minima:
    std::nth_element(active_cells.begin(),
                active_cells.begin() + upper_median_rank,
                active_cells.end(),
                boost::bind(
                  t.compare_strictly(),
                  boost::bind(Cell_min<Cell>(), _1),
		  boost::bind(Cell_min<Cell>(), _2)));
    
    Cell_iterator lower_median_cell =
      active_cells.begin() + upper_median_rank;
    Value lower_median = (lower_median_cell->min)();
    
    // compute lower median of cell's maxima:
    std::nth_element(active_cells.begin(),
                active_cells.begin() + lower_median_rank,
                active_cells.end(),
                boost::bind(
                  t.compare_strictly(),
                  boost::bind(Cell_max< Cell >(ccd), _1),
                  boost::bind(Cell_max< Cell >(ccd), _2)));
    
    Cell_iterator upper_median_cell =
      active_cells.begin() + lower_median_rank;
    Value upper_median = (upper_median_cell->max)(ccd);
    
    // restore lower_median_cell, if it has been displaced
    // by the second search
    if ((lower_median_cell->min)() != lower_median)
      lower_median_cell =
        find_if(active_cells.begin(),
                active_cells.end(),
                boost::bind(
		  equal_to< Value >(), 
		  lower_median,
		  boost::bind(Cell_min< Cell >(), _1)));
    CGAL_optimisation_assertion(lower_median_cell != active_cells.end());
    // ------------------------------------------------------
    // test feasibility of medians and remove cells accordingly:
    Cell_iterator new_end;
    
    
    if ( t.is_feasible( lower_median))
      if ( t.is_feasible( upper_median)) {
        // lower_median and upper_median are feasible
    
        // discard cells with all entries greater than
        // min( lower_median, upper_median) except for
        // one cell defining this minimum
    
        Cell_iterator min_median_cell;
        Value min_median;
        if ( lower_median < upper_median) {
          min_median_cell = lower_median_cell;
          min_median = lower_median;
        }
        else {
          min_median_cell = upper_median_cell;
          min_median = upper_median;
        }
    
        // save min_median_cell:
        iter_swap( min_median_cell, active_cells.begin());
    
        new_end =
          remove_if(
            active_cells.begin() + 1,
            active_cells.end(),
            boost::bind(
              t.compare_non_strictly(), 
	      min_median,
              boost::bind(Cell_min< Cell >(), _1)));
    
      } // lower_median and upper_median are feasible
      else { // lower_median is feasible, but upper_median is not
    
        // discard cells with all entries greater than
        // lower_median or all entries smaller than
        // upper_median except for the lower median cell
    
        // save lower_median_cell:
        iter_swap( lower_median_cell, active_cells.begin());
    
        new_end =
          remove_if(
            active_cells.begin() + 1,
            active_cells.end(),
            boost::bind(
              logical_or< bool >(),
	      boost::bind(
	        t.compare_non_strictly(),
		lower_median,
                boost::bind(Cell_min< Cell >(), _1)),
	      boost::bind(
                  t.compare_non_strictly(),
                  boost::bind(Cell_max< Cell >( ccd), _1),
		  upper_median)));
    
      } // lower_median is feasible, but upper_median is not
    else
      if ( t.is_feasible( upper_median)) {
        // upper_median is feasible, but lower_median is not
    
        // discard cells with all entries greater than
        // upper_median or all entries smaller than
        // lower_median except for the upper median cell
    
        // save upper_median_cell:
        iter_swap( upper_median_cell, active_cells.begin());
    
        new_end =
          remove_if(
            active_cells.begin() + 1,
            active_cells.end(),
	    boost::bind(
              logical_or< bool >(),
	      boost::bind(
		t.compare_non_strictly(),
		upper_median,
                boost::bind(Cell_min< Cell >(), _1)),
	      boost::bind(
                t.compare_non_strictly(),
		boost::bind(Cell_max< Cell >( ccd), _1),
		lower_median)));
    
      } // upper_median is feasible, but lower_median is not
      else { // both upper_median and lower_median are infeasible
    
        // discard cells with all entries smaller than
        // max BOOST_PREVENT_MACRO_SUBSTITUTION ( lower_median, upper_median)
    
        new_end =
          remove_if(
            active_cells.begin(),
            active_cells.end(),
            boost::bind(
              t.compare_non_strictly(),
	      boost::bind(Cell_max< Cell >( ccd), _1),
	      max BOOST_PREVENT_MACRO_SUBSTITUTION ( lower_median, upper_median)));
    
      } // both upper_median and lower_median are infeasible
    
      active_cells.erase( new_end, active_cells.end());
  } // for (;;)

  // there must be only one cell left:
  CGAL_optimisation_assertion( active_cells.size() == 1);
  CGAL_optimisation_assertion( ccd == 1);

  return ((*active_cells.begin()).min)();
}

} //namespace CGAL

#endif // ! (CGAL_SORTED_MATRIX_SEARCH_H)
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------
