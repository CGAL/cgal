// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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

#ifndef CGAL_POLYNOMIAL_INTERNAL_DERIVATIVE_H
#define CGAL_POLYNOMIAL_INTERNAL_DERIVATIVE_H

#include <CGAL/Polynomial/basic.h>
#include <iterator>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class Fn>
struct Derivative
{
private:
  typedef typename Fn::iterator Cit;
  struct It
  {
    It(Cit it, int i): i_(i), cit_(it) {}
    It(){}
    typedef typename Fn::iterator MCit;
    typedef typename std::iterator_traits<MCit> Traits;
    typedef typename Traits::iterator_category iterator_category;
    typedef typename Traits::value_type value_type;
    typedef typename Traits::pointer pointer;
    typedef typename Traits::reference reference;
    typedef typename Traits::difference_type difference_type;

    value_type operator*() const
    {
      return value_type(i_)**cit_;
    }
    pointer operator->() const
    {
      return value_type(i_)**cit_;
    }	  
    bool operator<(It o) const
    {
      return cit_ < o.cit_;
    }
    It operator++() {
      ++i_;
      ++cit_;
      return *this;
    }
    It operator++(int) {
      It t=*this;
      ++i_;
      ++cit_;
      return t;
    }
    It operator--() {
      --i_;
      --cit_;
      return *this;
    }
    It operator--(int) {
      It t=*this;
      --i_;
      --cit_;
      return t;
    }
    bool operator==(It o) const
    {
      return cit_ == o.cit_;
    }
    bool operator!=(It o) const
    {
      return cit_ != o.cit_;
    }
    difference_type operator-(It o) const
    {
      return cit_-o.cit_;
    }
    It operator+=(difference_type i) {
      cit_+= i;
      i_+= i;
      return *this;
    }
  protected:
    int i_;
    MCit cit_;
  };
public:
  typedef Fn result_type;
  typedef Fn argument_type;

  result_type operator()(const argument_type &o) const
  {
    if (o.is_constant()) { return result_type(typename result_type::NT(0)); }
    else {
      It b(o.begin(),0);
      ++b;
      return result_type(b, It(o.end(), o.degree()+1));
    }
  }
  
  void write(std::ostream &out) const
  {
    out << "diff";
  }
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
