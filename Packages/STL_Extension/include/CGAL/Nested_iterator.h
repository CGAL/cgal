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
// file          : include/CGAL/Nested_iterator.h
// package       : STL_Extension
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================


#ifndef CGAL_NESTED_ITERATOR_H
#define CGAL_NESTED_ITERATOR_H

#include <CGAL/iterator.h>

CGAL_BEGIN_NAMESPACE

template<class It>
struct Nested_iterator_traits
{
  typedef It                                 Base_iterator;
  typedef typename It::value_type::iterator  Nested_iterator;

  Nested_iterator begin(It it) const { return it->begin(); }
  Nested_iterator end(It it) const { return it->end(); }
};

namespace CGALi {

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
    typedef typename Predicate::Traits::Nested_iterator  Nested_iterator;
    typedef typename Predicate::Traits::Base_iterator    Base_iterator;

  public:    
    FI_w_begin_end() : F_iterator() {}

    FI_w_begin_end(Base_iterator it1,
		   Base_iterator it2,
		   Base_iterator it3)
      : F_iterator(it1, it2, Predicate(), it3) {}

    Nested_iterator begin(Base_iterator it)
    {
      return predicate().begin(it);
    }
    Nested_iterator end(Base_iterator it)
    {
      return predicate().end(it);
    }

  };

}

template <typename Base_it,
	  typename Tr = Nested_iterator_traits<Base_it> >
class Nested_iterator
  : private CGALi::FI_w_begin_end<
  Filter_iterator<Base_it, CGALi::Emptyness_predicate<Tr> > >
{
public:
  typedef Base_it                               Base_iterator;
  typedef Tr                                    Traits;
  typedef typename Tr::Nested_iterator          Iterator;

protected:
  typedef Nested_iterator<Base_iterator, Traits>  Self;
  typedef CGALi::
  FI_w_begin_end< Filter_iterator<Base_iterator,
				  CGALi::Emptyness_predicate<Tr> > >
  Filter_base_iterator;

public:
  typedef typename Iterator::reference          reference;
  typedef typename Iterator::pointer            pointer;
  typedef typename Iterator::value_type         value_type;
  typedef typename Iterator::difference_type    difference_type;
  typedef typename Iterator::iterator_category  iterator_category;

public:
  Nested_iterator() : Filter_base_iterator(), nested_it_() {}
  Nested_iterator(Base_iterator base_it_begin,
		  Base_iterator base_it_end,
		  Base_iterator base_it_cur)
    : Filter_base_iterator(base_it_begin, base_it_end, base_it_cur),
      nested_it_()
  {
    if ( !this->is_end() ) {
      nested_it_ = this->begin( this->base() );
    }
  }

  Self& operator++()
  {
    if ( nested_it_ != end( this->base() ) ) {
      ++nested_it_;
      if ( nested_it_ == end( this->base() ) ) {
	Filter_base_iterator::operator++();
	if ( !this->is_end() ) {
	  nested_it_ = begin( this->base() );
	}
      }
    } 
    return *this;
  }

  Self& operator++(int)
  {
    Self tmp = *this;
    ++(*this);
    return tmp;
  }

  Self& operator--()
  {
    if ( this->is_end() ) {
      Filter_base_iterator::operator--();
      nested_it_ = --end(this->base());
    } else {
      if ( nested_it_ != begin( this->base() ) ) {
	--nested_it_;
      } else {
	Filter_base_iterator::operator--();
	nested_it_ = --end( this->base() );
      }
    }
    return *this;
  }

  Self& operator--(int)
  {
    Self tmp = *this;
    --(*this);
    return tmp;
  }

  
  reference  operator*()  const
  {
    return nested_it_.operator*();
  }

  pointer    operator->() const
  {
    return nested_it_.operator->();
  }

  friend bool operator== <>(const Self&, const Self&);

protected:
  Iterator       nested_it_;
};





template<class Base_it, class Traits>
inline
bool operator==(const Nested_iterator<Base_it,Traits>& it1,
		const Nested_iterator<Base_it,Traits>& it2)
{
  CGAL_precondition( it1.b_ == it2.b_ && it1.e_ == it2.e_ );

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


CGAL_END_NAMESPACE

#endif // CGAL_NESTED_ITERATOR_H

	
