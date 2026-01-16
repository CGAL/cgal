// Copyright (c) 2003
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>


#ifndef CGAL_CONCATENATE_ITERATOR_H
#define CGAL_CONCATENATE_ITERATOR_H

#include <CGAL/basic.h>
#include <iterator>
#include <type_traits>


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
  typedef std::iterator_traits<It2>            Traits2;

public:
  typedef It1                                  Iterator1;
  typedef It2                                  Iterator2;

  typedef typename Traits1::reference          reference;
  typedef typename Traits1::pointer            pointer;
  typedef typename Traits1::value_type         value_type;
  typedef typename Traits1::difference_type    difference_type;
  typedef std::common_type_t<
    typename Traits1::iterator_category,
    typename Traits2::iterator_category
  > iterator_category;

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

  Self operator+(difference_type offset) const
  {
    if (offset < 0) {
      return this->operator-(-offset);
    }
    if constexpr (std::is_same_v<iterator_category, std::random_access_iterator_tag>) {
      // O(1) for random access iterators
      if (i1_ != e1_) {
        difference_type remaining_in_first = e1_ - i1_;
        if (offset < remaining_in_first) {
          return Self(e1_, b2_, i1_ + offset);
        } else {
          return Self(e1_, b2_, b2_ + (offset - remaining_in_first), 0);
        }
      } else {
        return Self(e1_, b2_, i2_ + offset, 0);
      }
    } else {
      // O(n) fallback for non-random-access iterators
      Self res(*this);
      for (difference_type i = 0; i < offset; ++i) { ++res; }
      return res;
    }
  }

  Self& operator+=(difference_type offset)
  {
    *this = this->operator+(offset);
    return *this;
  }

  Self operator-(difference_type offset) const
  {
    if (offset < 0) {
      return this->operator+(-offset);
    }
    if constexpr (std::is_same_v<iterator_category, std::random_access_iterator_tag>) {
      // O(1) for random access iterators
      if (i1_ != e1_) {
        // Currently in first range, just go back
        return Self(e1_, b2_, i1_ - offset);
      } else {
        difference_type pos_in_second = i2_ - b2_;
        if (offset <= pos_in_second) {
          return Self(e1_, b2_, i2_ - offset, 0);
        } else {
          // Need to go back into first range
          return Self(e1_, b2_, e1_ - (offset - pos_in_second));
        }
      }
    } else {
      // O(n) fallback for non-random-access iterators
      Self res(*this);
      for (difference_type i = 0; i < offset; ++i) { --res; }
      return res;
    }
  }

  Self& operator-=(difference_type offset)
  {
    *this = this->operator-(offset);
    return *this;
  }

  difference_type operator-(const Self& other) const
  {
    if constexpr (std::is_same_v<iterator_category, std::random_access_iterator_tag>) {
      // O(1) for random access iterators
      bool this_in_first = (i1_ != e1_);
      bool other_in_first = (other.i1_ != other.e1_);

      if (this_in_first && other_in_first) {
        return i1_ - other.i1_;
      } else if (!this_in_first && !other_in_first) {
        return i2_ - other.i2_;
      } else if (!this_in_first && other_in_first) {
        return (e1_ - other.i1_) + (i2_ - b2_);
      } else {
        return -((other.e1_ - i1_) + (other.i2_ - other.b2_));
      }
    } else {
      // operator-(iterator) is only valid for random access iterators
      // Return 0 as a fallback (undefined behavior for non-random-access)
      CGAL_assertion_msg(false, "operator-(iterator) requires random access iterators");
      return difference_type(0);
    }
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
