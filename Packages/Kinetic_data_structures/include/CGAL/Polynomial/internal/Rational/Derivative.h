// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_INTERNAL_DERIVATIVE_H
#define CGAL_POLYNOMIAL_INTERNAL_DERIVATIVE_H

#include <CGAL/Polynomial/basic.h>
#include <vector>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <class Fn>
struct Derivative{
private:
  typedef typename Fn::iterator Cit;
  struct It {
    It(Cit it, int i): i_(i), cit_(it) {}
    It(){}
    typedef typename Fn::iterator MCit;
    typedef typename MCit::iterator_category iterator_category;
    typedef typename MCit::value_type value_type;
    typedef typename MCit::pointer pointer;
    typedef typename MCit::reference reference;
    typedef typename MCit::difference_type difference_type;

    value_type operator*() const {
      return value_type(i_)**cit_;
    }
    pointer operator->() const {
      return value_type(i_)**cit_;
    }
    It operator++(){
      ++i_;
      ++cit_;
      return *this;
    }
    It operator++(int){
      It t=*this;
      ++i_;
      ++cit_;
      return t;
    }
    It operator--(){
      --i_;
      --cit_;
      return *this;
    }
    It operator--(int){
      It t=*this;
      --i_;
      --cit_;
      return t;
    }
    bool operator==(It o) const {
      return cit_ == o.cit_;
    }
    bool operator!=(It o) const {
      return cit_ != o.cit_;
    }
    difference_type operator-(It o) const {
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

  result_type operator()(const argument_type &o) const {
    if (o.is_constant()) { return result_type(typename result_type::NT(0)); }
    else {
      It b(o.begin(),0);
      ++b;
      return result_type(b, It(o.end(), o.degree()+1));
    }
  }

  void write(std::ostream &out) const {
    out << "diff";
  }
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
