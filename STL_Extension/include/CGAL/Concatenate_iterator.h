// Copyright (c) 2003  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>


#ifndef CGAL_CONCATENATE_ITERATOR_H
#define CGAL_CONCATENATE_ITERATOR_H

#include <CGAL/basic.h>
#include <iterator>


#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4396)
#endif
namespace CGAL {

template <class It1, class It2> class Concatenate_iterator;

template <class It1, class It2>
bool operator==(const Concatenate_iterator<It1,It2>&,
		const Concatenate_iterator<It1,It2>&);


template <class It1, class It2>
class Concatenate_iterator
{
private:
  typedef Concatenate_iterator<It1,It2>        Self;
  typedef std::iterator_traits<It1>            Traits1;

public:
  typedef It1                                  Iterator1;
  typedef It2                                  Iterator2;

  typedef typename Traits1::reference          reference;
  typedef typename Traits1::pointer            pointer;
  typedef typename Traits1::value_type         value_type;
  typedef typename Traits1::difference_type    difference_type;
  typedef typename Traits1::iterator_category  iterator_category;

public:
  Concatenate_iterator() : e1_(), i1_(), b2_(), i2_() {}
  Concatenate_iterator(It1 e1, It2 b2, It1 i1)
    : e1_(e1), i1_(i1), b2_(b2), i2_(b2) {}
  Concatenate_iterator(It1 e1, It2 b2, It2 i2, int)
    : e1_(e1), i1_(e1), b2_(b2), i2_(i2) {}

  Self& operator++()
  {
    if ( i1_ == e1_ ) {
      ++i2_;
    } else {
      ++i1_;
    }
    return *this;
  }

  Self operator++(int)
  {
    Self tmp = *this;
    ++(*this);
    return tmp;
  }

  Self& operator--()
  {
    if ( i2_ == b2_ ) {
      --i1_;
    } else {
      --i2_;
    }
    return *this;
  }

  Self operator--(int)
  {
    Self tmp = *this;
    --(*this);
    return tmp;
  }

  
  reference  operator*()  const
  {
    if ( i1_ == e1_ ) {
      return *i2_;
    } else {
      return *i1_;
    }
  }

  pointer    operator->() const
  {
    if ( i1_ == e1_ ) {
      return i2_.operator->();
    } else {
      return i1_.operator->();
    }
  }

  friend bool operator==<>(const Self&, const Self&);

protected:
  It1 e1_, i1_;
  It2 b2_, i2_;
};



template<class It1, class It2>
inline
bool operator==(const Concatenate_iterator<It1, It2>& it1,
		const Concatenate_iterator<It1, It2>& it2)
{
  return (it1.i1_ == it2.i1_ && it1.i2_ == it2.i2_);
}

template<class It1, class It2>
inline
bool operator!=(const Concatenate_iterator<It1, It2>& it1,
		const Concatenate_iterator<It1, It2>& it2)
{
  return !(it1 == it2);
}


} //namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif



#endif // CGAL_CONCATENATE_ITERATOR
