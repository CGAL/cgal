// Copyright (c) 2002  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_TRIVIAL_ITERATOR_H
#define CGAL_TRIVIAL_ITERATOR_H

#include <CGAL/license/Triangulation_2.h>


#include <iterator>
#include <CGAL/iterator.h>

namespace CGAL { 

// TODO :
// - comparison operators should be global, but it causes problems...
// - Have a look at Boost's concept_checking and archetypes :
//   http://www.boost.org/libs/concept_check/concept_check.htm

class Trivial_iterator_tag{};


template <class I>
class Trivial_iterator
{
public:
  typedef I                                                 Iterator;
  typedef Trivial_iterator<I>                               Self;
  typedef typename std::iterator_traits<I>::value_type      value_type;
  typedef typename std::iterator_traits<I>::difference_type difference_type;
  typedef typename std::iterator_traits<I>::reference       reference;
  typedef typename std::iterator_traits<I>::pointer         pointer;
  typedef Trivial_iterator_tag                              iterator_category;
  // Special for circulators.
  typedef I_Circulator_size_traits<iterator_category,I>     C_S_Traits;
  typedef typename  C_S_Traits::size_type                   size_type;


  Trivial_iterator() {}
  Trivial_iterator(const I &i) : base_(i) {}

  // To allow conversion from iterator to const_iterator.
  template <class Iter>
  Trivial_iterator(const Trivial_iterator<Iter> &t)
    : base_(t.base()) {}

  reference operator*() const { return *base_;  }
  pointer operator->() const  { return &*base_; }

  bool operator==(const Trivial_iterator &b) const { return base()==b.base(); }
  bool operator!=(const Trivial_iterator &b) const { return base()!=b.base(); }

private:
  const Iterator & base() const { return base_; }

  Iterator base_;
};



class Trivial_comparable_iterator_tag{};

template <class I>
class Trivial_comparable_iterator
{
public:
  typedef I                                            Iterator;
  typedef Trivial_comparable_iterator<I>               Self;
  typedef typename std::iterator_traits<I>::value_type value_type;
  typedef typename std::iterator_traits<I>::difference_type
                                                       difference_type;
  typedef typename std::iterator_traits<I>::reference  reference;
  typedef typename std::iterator_traits<I>::pointer    pointer;
  typedef Trivial_comparable_iterator_tag              iterator_category;
  // Special for circulators.
  typedef I_Circulator_size_traits<iterator_category,I> C_S_Traits;
  typedef typename  C_S_Traits::size_type               size_type;


  Trivial_comparable_iterator() {}
  Trivial_comparable_iterator(const I &i) : base_(i) {}

  // To allow conversion from iterator to const_iterator.
  template <class Iter>
  Trivial_comparable_iterator(const Trivial_comparable_iterator<Iter> &t)
    : base_(t.base()) {}

  reference operator*() const { return *base_;  }
  pointer operator->() const  { return &*base_; }

  bool operator==(const Trivial_comparable_iterator &b) const 
    { return base()==b.base(); }
  bool operator!=(const Trivial_comparable_iterator &b) const  
    { return base()!=b.base(); }

  bool operator< (const Trivial_comparable_iterator &b) const  
    { return base()< b.base(); }
  bool operator> (const Trivial_comparable_iterator &b) const  
    { return base()> b.base(); }
  bool operator<=(const Trivial_comparable_iterator &b) const  
    { return base()<=b.base(); }
  bool operator>=(const Trivial_comparable_iterator &b) const  
    { return base()>=b.base(); }

private:
  const Iterator & base() const { return base_; }

  Iterator base_;
};


// Some macros depending on CGAL_NO_CONCEPT_CHECKING.
#ifndef CGAL_NO_CONCEPT_CHECKING
#  define CGAL_TRIVIAL_ITERATOR_CHECKER(X)    CGAL::Trivial_iterator<X>
#  define CGAL_TRIVIAL_COMPARABLE_ITERATOR_CHECKER(X) \
     CGAL::Trivial_comparable_iterator<X>
#else
#  define CGAL_TRIVIAL_ITERATOR_CHECKER(X)    X
#  define CGAL_TRIVIAL_COMPARABLE_ITERATOR_CHECKER(X) X
#endif

} //namespace CGAL

#endif // CGAL_TRIVIAL_ITERATOR_H
