// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
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
// file          : include/CGAL/Trivial_iterator.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//
// ======================================================================

#ifndef CGAL_TRIVIAL_ITERATOR_H
#define CGAL_TRIVIAL_ITERATOR_H

#include <iterator>
#include <CGAL/iterator.h>

CGAL_BEGIN_NAMESPACE 

// TODO :
// - comparison operators should be global, but it causes problems...
// - Have a look at Boost's concept_checking and archetypes :
//   http://www.boost.org/libs/concept_check/concept_check.htm

class Trivial_iterator_tag{};

#if defined CGAL_CFG_NO_ITERATOR_TRAITS && \
   !defined CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
template < class I, class Ref, class Ptr, class Val, class Dist >
#else
template < class I>
#endif
class Trivial_iterator
{
public:
  typedef I                                                 Iterator;
#if defined  CGAL_CFG_NO_ITERATOR_TRAITS || defined _MSC_VER
  typedef Trivial_iterator<I,Ref,Ptr,Val,Dist>              Self;
  typedef Val                                               value_type;
  typedef Dist                                              difference_type;
  typedef Ref                                               reference;
  typedef Ptr                                               pointer;
#else
  typedef Trivial_iterator<I>                               Self;
  typedef typename std::iterator_traits<I>::value_type      value_type;
  typedef typename std::iterator_traits<I>::difference_type difference_type;
  typedef typename std::iterator_traits<I>::reference       reference;
  typedef typename std::iterator_traits<I>::pointer         pointer;
#endif
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

#if defined CGAL_CFG_NO_ITERATOR_TRAITS && \
   !defined CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
template < class I, class Ref, class Ptr, class Val, class Dist >
#else
template < class I>
#endif
class Trivial_comparable_iterator
{
public:
  typedef I                                            Iterator;
#if defined  CGAL_CFG_NO_ITERATOR_TRAITS || defined _MSC_VER
  typedef Trivial_comparable_iterator<I,Ref,Ptr,Val,Dist>  Self;
  typedef Val                                          value_type;
  typedef Dist                                         difference_type;
  typedef Ref                                          reference;
  typedef Ptr                                          pointer;
#else
  typedef Trivial_comparable_iterator<I>               Self;
  typedef typename std::iterator_traits<I>::value_type value_type;
  typedef typename std::iterator_traits<I>::difference_type
                                                       difference_type;
  typedef typename std::iterator_traits<I>::reference  reference;
  typedef typename std::iterator_traits<I>::pointer    pointer;
#endif
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

// This macro added to workaround the fact that MSC does not support
// default template parameter
#if defined _MSC_VER && !defined CGAL_NO_CONCEPT_CHECKING
#  define CGAL_TRIVIAL_ITERATOR_CHECKER_POINTER(X) \
          CGAL::Trivial_iterator<X*, X&, X*, X, std::ptrdiff_t>
#  define CGAL_TRIVIAL_COMPARABLE_ITERATOR_CHECKER_POINTER(X) \
          CGAL::Trivial_comparable_iterator<X*, X&, X*, X, std::ptrdiff_t>
#else
#  define CGAL_TRIVIAL_ITERATOR_CHECKER_POINTER(X) \
          CGAL_TRIVIAL_ITERATOR_CHECKER(X*)
#  define CGAL_TRIVIAL_COMPARABLE_ITERATOR_CHECKER_POINTER(X) \
          CGAL_TRIVIAL_COMPARABLE_ITERATOR_CHECKER(X*)
#endif


// For backward compatibility (can be removed soon, it's not been release)
#define CGAL_COMPARABLE_ITERATOR_CHECKER(X) \
        CGAL_TRIVIAL_COMPARABLE_ITERATOR_CHECKER(X)
#define CGAL_COMPARABLE_ITERATOR_CHECKER_POINTER(X) \
        CGAL_TRIVIAL_COMPARABLE_ITERATOR_CHECKER_POINTER(X)

CGAL_END_NAMESPACE

#endif // CGAL_TRIVIAL_ITERATOR_H
