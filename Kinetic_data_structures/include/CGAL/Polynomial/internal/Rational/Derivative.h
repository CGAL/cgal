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

#ifndef CGAL_POLYNOMIAL_INTERNAL_DERIVATIVE_H
#define CGAL_POLYNOMIAL_INTERNAL_DERIVATIVE_H

#include <CGAL/Polynomial/basic.h>
#include <iterator>
#include <boost/iterator/iterator_facade.hpp>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class Fn>
struct Derivative
{
private:
  typedef typename Fn::iterator Cit;
  struct It :     
    public boost::iterator_facade<It, 
                                  typename std::iterator_traits<typename Fn::iterator>::value_type,
                                  typename std::iterator_traits<typename Fn::iterator>::iterator_category,
                                  // Use value_type as reference (see
                                  // dereference() below()) 
                                  typename std::iterator_traits<typename Fn::iterator>::value_type,
                                  typename std::iterator_traits<typename Fn::iterator>::difference_type>
  {
    It(Cit it, int i): i_(i), cit_(it) {}
    It(){}

  private:
    typedef typename Fn::iterator MCit;
    typedef typename std::iterator_traits<MCit> Traits;

  public:
    typedef typename Traits::difference_type difference_type;
    typedef typename Traits::value_type value_type;

  private:
    friend class boost::iterator_core_access;
    value_type dereference() const { return value_type(i_)**cit_; }
    void increment() { ++i_; ++cit_;}
    bool equal(const It& i) const { return (i_ == i.i_) && (cit_ == i.cit_); }
    void decrement() { --i_; --cit_;}
    void advance(const difference_type n) { i_ += n; cit_ += n; }
    difference_type distance_to(const It& o) const { return o.cit_ - cit_; }

  protected:
    int i_;
    MCit cit_;
  }; // end nested struct It

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
