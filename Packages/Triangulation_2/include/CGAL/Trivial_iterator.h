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
// - iterator_category should not be exactly forwarded...
// - comparison operators should be global.
// - Have a look at Boost's concept_checking and archetypes :
//   http://www.boost.org/libs/concept_check/concept_check.htm

#if defined(CGAL_CFG_NO_ITERATOR_TRAITS) && \
!defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
template < class I, class Ref, class Ptr,
           class Val, class Dist, class Ctg >
#elif !defined __SUNPRO_CC
template < class I,
           class Ref  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::reference,
           class Ptr  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::pointer,
           class Val  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::value_type,
           class Dist = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::difference_type,
           class Ctg  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::iterator_category>
#else
template < class I,
           class Ref  = typename CGALi::IT_rename<I>::REF,
           class Ptr  = typename CGALi::IT_rename<I>::PTR,
           class Val  = typename CGALi::IT_rename<I>::VAL,
           class Dist = typename CGALi::IT_rename<I>::DIF,
           class Ctg  = typename CGALi::IT_rename<I>::CAT >
#endif
class Trivial_iterator
{
public:
  typedef I                                            Iterator;
  typedef Trivial_iterator<I,Ref,Ptr,Val,Dist,Ctg>     Self;
  typedef Ctg                                          iterator_category;
  typedef Val                                          value_type;
  typedef Dist                                         difference_type;
#ifdef CGAL_CFG_NO_ITERATOR_TRAITS
  typedef Ref                                          reference;
  typedef Ptr                                          pointer;
#else
  typedef typename std::iterator_traits<I>::reference  reference;
  typedef typename std::iterator_traits<I>::pointer    pointer;
#endif
  // Special for circulators.
  typedef I_Circulator_size_traits<iterator_category,I> C_S_Traits;
  typedef typename  C_S_Traits::size_type               size_type;


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

#if defined(CGAL_CFG_NO_ITERATOR_TRAITS) && \
!defined(CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT)
template < class I, class Ref, class Ptr,
           class Val, class Dist, class Ctg >
#elif !defined __SUNPRO_CC
template < class I,
           class Ref  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::reference,
           class Ptr  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::pointer,
           class Val  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::value_type,
           class Dist = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::difference_type,
           class Ctg  = CGAL_TYPENAME_MSVC_NULL
                        std::iterator_traits<I>::iterator_category>
#else
template < class I,
           class Ref  = typename CGALi::IT_rename<I>::REF,
           class Ptr  = typename CGALi::IT_rename<I>::PTR,
           class Val  = typename CGALi::IT_rename<I>::VAL,
           class Dist = typename CGALi::IT_rename<I>::DIF,
           class Ctg  = typename CGALi::IT_rename<I>::CAT >
#endif
class Comparable_iterator
{
public:
  typedef I                                            Iterator;
  typedef Comparable_iterator<I,Ref,Ptr,Val,Dist,Ctg>     Self;
  typedef Ctg                                          iterator_category;
  typedef Val                                          value_type;
  typedef Dist                                         difference_type;
#if defined  CGAL_CFG_NO_ITERATOR_TRAITS || defined _MSC_VER
  typedef Ref                                          reference;
  typedef Ptr                                          pointer;
#else
  typedef typename std::iterator_traits<I>::reference  reference;
  typedef typename std::iterator_traits<I>::pointer    pointer;
#endif
  // Special for circulators.
  typedef I_Circulator_size_traits<iterator_category,I> C_S_Traits;
  typedef typename  C_S_Traits::size_type               size_type;


  Comparable_iterator() {}
  Comparable_iterator(const I &i) : base_(i) {}

  // To allow conversion from iterator to const_iterator.
  template <class Iter>
  Comparable_iterator(const Comparable_iterator<Iter> &t)
    : base_(t.base()) {}

  reference operator*() const { return *base_;  }
  pointer operator->() const  { return &*base_; }

  bool operator==(const Comparable_iterator &b) const 
    { return base()==b.base(); }
  bool operator!=(const Comparable_iterator &b) const  
    { return base()!=b.base(); }

  bool operator< (const Comparable_iterator &b) const  
    { return base()< b.base(); }
  bool operator> (const Comparable_iterator &b) const  
    { return base()> b.base(); }
  bool operator<=(const Comparable_iterator &b) const  
    { return base()<=b.base(); }
  bool operator>=(const Comparable_iterator &b) const  
    { return base()>=b.base(); }

private:
  const Iterator & base() const { return base_; }

  Iterator base_;
};

class Comparable_iterator_tag{};

// This macro added to workaround the fact that MSC does not support
// defaut  template parameter
#ifdef _MSC_VER
#  ifdef CGAL_NO_CONCEPT_CHECKING 
#    define CGAL_COMPARABLE_ITERATOR_CHECKER_POINTER(X) X*
#  else
#    define CGAL_COMPARABLE_ITERATOR_CHECKER_POINTER(X) \
            CGAL::Comparable_iterator<X*, X&, X*, X, std::ptrdiff_t, \
                                      CGAL::Comparable_iterator_tag>
#  endif
#else
#  define CGAL_COMPARABLE_ITERATOR_CHECKER_POINTER(X) \
          CGAL_COMPARABLE_ITERATOR_CHECKER(X*)
#endif

#ifndef CGAL_NO_CONCEPT_CHECKING
#  define CGAL_TRIVIAL_ITERATOR_CHECKER(X)    CGAL::Trivial_iterator<X>
#  define CGAL_COMPARABLE_ITERATOR_CHECKER(X) CGAL::Comparable_iterator<X>
#else
#  define CGAL_TRIVIAL_ITERATOR_CHECKER(X)    X
#  define CGAL_COMPARABLE_ITERATOR_CHECKER(X) X
#endif

CGAL_END_NAMESPACE

#endif // CGAL_TRIVIAL_ITERATOR_H
