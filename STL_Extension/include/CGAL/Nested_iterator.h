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


#ifndef CGAL_NESTED_ITERATOR_H
#define CGAL_NESTED_ITERATOR_H

#include <CGAL/iterator.h>
#include <iterator>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4396)
#endif

namespace CGAL {

template<class It>
struct Nested_iterator_traits
{
  typedef It                                 Base_iterator;
  typedef typename std::iterator_traits<It>::value_type::iterator  Iterator;

  Iterator begin(It it) const { return it->begin(); }
  Iterator end(It it) const { return it->end(); }
};

namespace internal {

  template<class Tr>
  struct Emptyness_predicate : public Tr
  {
    typedef Tr            Traits;

    bool operator()(typename Traits::Base_iterator base_it) const
    {
      return this->begin(base_it) == this->end(base_it);
    }


  };

  template<class F_iterator>
  class FI_w_begin_end : public F_iterator
  {
  private:
    typedef typename F_iterator::Predicate               Predicate;
    typedef typename Predicate::Traits::Iterator         Iterator;
    typedef typename Predicate::Traits::Base_iterator    Base_iterator;

  public:
    FI_w_begin_end() : F_iterator() {}

    FI_w_begin_end(Base_iterator it2,
                   Base_iterator it3)
      : F_iterator(it2, Predicate(), it3) {}

    Iterator begin(Base_iterator it)
    {
      return this->predicate().begin(it);
    }

    Iterator end(Base_iterator it)
    {
      return this->predicate().end(it);
    }

  };

}

template<class Base_it, class Tr> class Nested_iterator;

template<class Base_it, class Tr>
bool operator==(const Nested_iterator<Base_it,Tr>&,
                const Nested_iterator<Base_it,Tr>&);

template <typename Base_it,
          typename Tr = Nested_iterator_traits<Base_it> >
class Nested_iterator
  : private internal::FI_w_begin_end<
  Filter_iterator<Base_it, internal::Emptyness_predicate<Tr> > >
{
  typedef Nested_iterator<Base_it,Tr>           Self;
public:
  typedef Base_it                               Base_iterator;
  typedef Tr                                    Traits;
  typedef typename Tr::Iterator                 Iterator;

protected:
  typedef internal::
  FI_w_begin_end< Filter_iterator<Base_iterator,
                                  internal::Emptyness_predicate<Tr> > >
  Filter_base_iterator;

private:
  typedef std::iterator_traits<Iterator>        ItTraits;

public:
  typedef typename ItTraits::reference          reference;
  typedef typename ItTraits::pointer            pointer;
  typedef typename ItTraits::value_type         value_type;
  typedef typename ItTraits::difference_type    difference_type;
  typedef typename ItTraits::iterator_category  iterator_category;

public:
  Nested_iterator() : Filter_base_iterator(), nested_it_() {}
  Nested_iterator(Base_iterator base_it_end,
                  Base_iterator base_it_cur)
    : Filter_base_iterator(base_it_end, base_it_cur), nested_it_()
  {
    if ( !this->is_end() ) {
      nested_it_ = this->begin( this->base() );
    }
  }

  Nested_iterator(const Self& other)
    : Filter_base_iterator(other)
  {
    if ( !other.is_end() )
      nested_it_ = other.nested_it_;
  }

  Self& operator=(const Self& other)
  {
    copy_from(other);
    return *this;
  }

  Self& operator++()
  {
    if ( nested_it_ != this->end( this->base() ) ) {
      ++nested_it_;
      if ( nested_it_ == this->end( this->base() ) ) {
        Filter_base_iterator::operator++();
        if ( !this->is_end() ) {
          nested_it_ = this->begin( this->base() );
        }
      }
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
    if ( this->is_end() ) {
      Filter_base_iterator::operator--();
      nested_it_ = this->end(this->base());
      --nested_it_;
    } else {
      if ( nested_it_ != this->begin( this->base() ) ) {
        --nested_it_;
      } else {
        Filter_base_iterator::operator--();
        nested_it_ = this->end( this->base() );
        --nested_it_;
      }
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
    return *nested_it_;
  }

  pointer    operator->() const
  {
    return nested_it_.operator->();
  }

  friend bool operator==<>(const Self&, const Self&);

protected:
  void copy_from(const Self& other)
  {
    Filter_base_iterator::operator=(other);
    if ( !other.is_end() ) {
      nested_it_ = other.nested_it_;
    }
  }

protected:
  Iterator       nested_it_;
};





template<class Base_it, class Traits>
inline
bool operator==(const Nested_iterator<Base_it,Traits>& it1,
                const Nested_iterator<Base_it,Traits>& it2)
{
  //  CGAL_precondition( it1.b_ == it2.b_ && it1.e_ == it2.e_ );

  if ( it1.base() != it2.base() ) { return false; }
  return it1.is_end() || ( it1.nested_it_ == it2.nested_it_ );
}

template<class Base_it, class Traits>
inline
bool operator!=(const Nested_iterator<Base_it,Traits>& it1,
                const Nested_iterator<Base_it,Traits>& it2)
{
  return !(it1 == it2);
}


} //namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_NESTED_ITERATOR_H
