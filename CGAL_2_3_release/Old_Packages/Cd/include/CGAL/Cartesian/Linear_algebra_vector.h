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
// file          : include/CGAL/Cartesian/Linear_algebra_vector.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_LINEAR_ALGEBRA_VECTOR_D_H
#define CGAL_CARTESIAN_LINEAR_ALGEBRA_VECTOR_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/d_tuple.h>
#include <algorithm>
#include <functional>
#include <iostream>

CGAL_BEGIN_NAMESPACE

template < class LA_ >
class LA_vectorCd;

template < class LA >
std::ostream &
operator<<(std::ostream &, const LA_vectorCd<LA> &);

template < class LA_ >
class LA_vectorCd : public Handle
{
  // The last thing the world needs is a Vector class. Nevertheless,
  // we have to provide one because no *generic* linear algebra toolkit
  // is available -- *sigh*
public:
  typedef LA_                              LA;
  typedef typename LA::FT                  FT;
  typedef typename LA::RT                  RT;
  typedef LA_vectorCd<LA>                  Self;
  typedef FT*                              iterator;
  typedef const FT*                        const_iterator;
  
  LA_vectorCd(int dim = 0);
  LA_vectorCd(const FT &a);
  LA_vectorCd(const FT &a, const FT &b);
  LA_vectorCd(const FT &a, const FT &b, const FT &c);
  LA_vectorCd(const FT &a, const FT &b, const FT &c, const FT &d);
  LA_vectorCd(const Handle &v);
  LA_vectorCd(const Self &v);
  template < class InputIterator >
  LA_vectorCd(InputIterator first, InputIterator last)
  {
    new_rep(last-first);
    std::copy(first,last,begin());
  }
  ~LA_vectorCd();

  Self&          operator=(const Self &v);

  bool           operator==(const Self &p) const;
  bool           operator!=(const Self &p) const;

  // Unary operators
  Self           operator+() const { return *this; }
  Self           operator-() const;
  // Binary operators
  Self           operator+(const Self &w) const;
  Self           operator-(const Self &w) const;
  FT             operator*(const Self &w) const;
  Self           operator*(const FT &c) const;
  Self           operator/(const FT &c) const;

  // Component access
  long           id()              const { return (long) ptr(); }
  int            dimension()       const { return ptr()->d; }

  const_iterator begin()           const { return ptr()->e; }
  const_iterator end()             const { return ptr()->e + dimension(); }
  const FT      &operator[](int i) const { return *(begin()+i); }

// protected:
  iterator       begin()                 { return ptr()->e; }
  iterator       end()                   { return ptr()->e + dimension(); }
  FT            &operator[](int i)       { return *(begin()+i); }

// debug:
  void  print() const { std::cout << *this; }

private:
  void new_rep(int d)              { PTR = new _d_tuple<FT>(d); }
  const _d_tuple<FT>* ptr()  const { return (const _d_tuple<FT>*)PTR; }
  _d_tuple<FT>*       ptr()        { return (_d_tuple<FT>*)PTR; }
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Linear_algebra_vector_d.C>
#endif 

#endif // CGAL_CARTESIAN_VECTOR_D_H
