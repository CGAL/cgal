// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_CONVERTERS_H
#define CGAL_POLYNOMIAL_CONVERTERS_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/nt_converters.h>
#include <iterator>

namespace CGAL { namespace POLYNOMIAL {

//! from P1 to P2
template<class P1, class P2, class Converter_t =
	 CGAL_POLYNOMIAL_NS::NT_converter<typename P1::NT,typename P2::NT> >
class Polynomial_converter
{
  class Iterator
  {
  public:
    Iterator(typename P1::iterator it, const Converter_t& fun): fun_(fun), it_(it){}
    Iterator(){}
    typedef typename std::iterator_traits<typename P1::iterator> Traits;
    typedef typename Traits::iterator_category iterator_category;
    typedef typename Traits::difference_type difference_type;
    typedef typename Converter_t::result_type value_type;
    typedef value_type reference;
    typedef const value_type *pointer;

    value_type operator*() const
    {
      return fun_(*it_);
    }
    /*pointer_type operator->() const {
      return fun_(P1::iterator::operator*());
      }*/
    difference_type operator-(Iterator o) const
    {
      return it_-o.it_;
    }

    bool operator<(const Iterator& o) const
    {
      return   ( it_-o.it_ ) <  0;
    }


    bool operator==(const Iterator &o) const
    {
      return it_==o.it_;
    }

    bool operator!=(const Iterator &o) const
    {
      return it_!=o.it_;
    }

    Iterator operator++() {
      ++it_;
      return *this;
    }


    Iterator operator++(int) {
      Iterator t= *this;
      ++it_;
      return t;
    }

    Iterator operator+=(difference_type i) {
      it_+=i;
      return *this;
    }

  protected:
    Converter_t fun_;
    typename P1::iterator it_;
  };

public:
  typedef P1           argument_type;
  typedef P2           result_type;
  typedef Converter_t  NT_converter;

  explicit Polynomial_converter(NT_converter cv = NT_converter()) : cv_(cv) {}

  P2 operator()(const P1& p1) const
  {
    Iterator b(p1.begin(), cv_);
    Iterator e(p1.end(), cv_);
    return P2(b,e);
  }
  /*typename P2::NT operator()(const typename P1::NT &p) const {
    return cv_(p);
    }*/
  const NT_converter& nt_converter() const
  {
    return cv_;
  }
protected:
  NT_converter cv_;
};

template <class F>
class Polynomial_identity_converter
{
public:
  Polynomial_identity_converter(){}
  typedef F argument_type;
  typedef F result_type;
  typedef Identity_converter<typename argument_type::NT> NT_converter;
  const F& operator()(const F&f) const
  {
    return f;
  }
  NT_converter nt_converter() const
  {
    return NT_converter();
  }
};

} } //namespace CGAL::POLYNOMIAL
#endif                                            // CGAL_POLYNOMIAL_CONVERTER_H
