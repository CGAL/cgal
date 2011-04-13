// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Linear_algebra_matrix.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_LINEAR_ALGEBRA_MATRIX_D_H
#define CGAL_CARTESIAN_LINEAR_ALGEBRA_MATRIX_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/d_tuple.h>
#include <algorithm>
#include <iterator>
#include <functional>
#include <utility>
#include <iostream>

CGAL_BEGIN_NAMESPACE

template < class LA_ >
class LA_matrixCd;

template < class LA >
std::ostream &
operator<<(std::ostream &, const LA_matrixCd<LA> &);

template < class LA_ >
class LA_matrixCd : public Handle
{
  // Matrix class
  // Internal representation is a single n*m vector of FTs
  // Rows are stored consecutively
public:
  typedef LA_                              LA;
  typedef typename LA::FT                  FT;
  typedef typename LA::RT                  RT;
  typedef typename LA::Vector              Vector;
  typedef LA_matrixCd<LA>                  Self;
  typedef FT*                              iterator;
  typedef FT*                              const_iterator;
  
  // Uninitialized constructors
  LA_matrixCd(int n = 0);
  LA_matrixCd(int n, int m);
  LA_matrixCd(const std::pair<int,int> &d);
  // Copy constructor
  LA_matrixCd(const Self &v);
  // Constructor from a set of row vectors
  template < class InputIterator >
  LA_matrixCd(int n, const InputIterator &first, const InputIterator &last)
    { // InputIterator::value_type has dimension(), begin() and end()
      // and begin() is itself an operator, with an operator*()
      // which returns a type that is assignable to a FT
      CGAL_kernel_precondition( last-first == n);
      int m = first->dimension();
      new_rep(n,m);
      InputIterator i;
      ptrdiff_t j;
      for (i=first,j=0; i!=last; ++i, j+=m)
        std::copy_n(i->begin(),m,begin()+j);
    }
  // Constructor from a single (long) vector of FTs
  template < class InputIterator >
  LA_matrixCd(int n, int m,
              const InputIterator &first, const InputIterator &last)
    { // InputIterator::value_type is assignable to a FT
      CGAL_kernel_precondition( last-first == n*m);
      new_rep(n,m);
      std::copy(first,last,begin());
    }
  ~LA_matrixCd();

  Self&          operator=(const Self &v);

  bool           operator==(const Self &p) const;
  bool           operator!=(const Self &p) const;

  // Unary operators
  Self           operator+() const { return *this; }
  Self           operator-() const;

  // Binary operators
  Self           operator+(const Self &w) const;
  Self           operator-(const Self &w) const;
  Self           operator*(const FT &c) const;
  Self           operator/(const FT &c) const;

  // Matrix and Matrix/vector multiplications
  Self           operator*(const Self &w) const;
  Vector         operator*(const Vector &w) const;

  std::pair<int,int> dimension() const; // (row,column)

  // Row and column
  int            row_dimension()    const { return _rdim; } // number of lines
  int            column_dimension() const { return _cdim; } // number of columns

  // Row and column copy access
  Vector         row(int i) const;
  Vector         column(int j) const;

  // Internal access
  const_iterator begin()             const { return ptr()->e; }
  const_iterator end()               const { return ptr()->e + ptr()->d; }
  
  const_iterator row_begin(int i)    const { return begin() + i*_cdim; }
  const_iterator row_end(int i)      const { return begin() + (i+1)*_cdim; }
  const_iterator operator[](int i)   const { return row_begin(i); }

  // we could return a special iterator type in which ++ is overloaded
  // to take into account the step by column_dimension()
  // if so, we could use std::copy() etc.
  // it's amazing there is no such adaptor for iterators
  // so we must step manually by column_dimension()
  const_iterator column_begin(int j) const { return begin() + j; }
  const_iterator column_end(int j)   const { return end() + j; }

  void swap_columns(int j, int k);

// protected:
  iterator       begin()                   { return ptr()->e; }
  iterator       end()                     { return ptr()->e + ptr()->d; }
  iterator       row_begin()               { return ptr()->e + i*_cdim; }
  iterator       row_end()                 { return ptr()->e + (i+1)*_cdim; }
  iterator       column_begin(int j)       { return begin() + j; }
  iterator       column_end(int j)         { return end() + j; }

// debug:
  void print() const { std::cout << *this; }

private:
  int _rdim, _cdim;
  void new_rep(int n, int m)               { _rdim = n; _cdim = m;
                                             PTR = new _d_tuple<FT>(n*m); }
  const _d_tuple<FT>* ptr()          const { return (const _d_tuple<FT>*)PTR; }
  _d_tuple<FT>*       ptr()                { return (_d_tuple<FT>*)PTR; }
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Linear_algebra_matrix_d.C>
#endif 

#endif // CGAL_CARTESIAN_MATRIX_D_H
