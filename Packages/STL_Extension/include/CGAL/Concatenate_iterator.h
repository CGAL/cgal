// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : include/CGAL/Concatenate_iterator.h
// package       : STL_Extension
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================


#ifndef CGAL_CONCATENATE_ITERATOR_H
#define CGAL_CONCATENATE_ITERATOR_H

#include <CGAL/basic.h>


CGAL_BEGIN_NAMESPACE


template <class It1, class It2>
class Concatenate_iterator
{
private:
  typedef Concatenate_iterator<It1,It2>    Self;

public:
  typedef It1                              Iterator1;
  typedef It2                              Iterator2;

  typedef typename It1::reference          reference;
  typedef typename It1::pointer            pointer;
  typedef typename It1::value_type         value_type;
  typedef typename It1::difference_type    difference_type;
  typedef typename It1::iterator_category  iterator_category;

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

  friend bool operator== <>(const Self&, const Self&);

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


CGAL_END_NAMESPACE




#endif // CGAL_CONCATENATE_ITERATOR

