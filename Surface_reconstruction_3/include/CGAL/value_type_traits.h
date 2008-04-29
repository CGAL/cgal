//=============================================================================
//
//  Traits class to get the value type of an output iterator.
//  Based on code posted by Alberto Ganesh Barbati at
//  http://www.adras.com/Why-no-std-back-insert-iterator-value-type.t2639-153-3.html
//
//=============================================================================

#ifndef CGAL_VALUE_TYPE_TRAITS_H
#define CGAL_VALUE_TYPE_TRAITS_H

#include <CGAL/basic.h>

#include <iterator>

CGAL_BEGIN_NAMESPACE


// Traits class to get the value type of any iterator,
// including an output iterator.
//
// Usage is:
// typedef typename value_type_traits<Iter>::type value_type;

template <class T>
struct value_type_traits
{
  typedef typename std::iterator_traits<T>::value_type type;
};

template <class T>
struct value_type_traits<std::back_insert_iterator<T> >
{
  typedef typename T::value_type type;
};


CGAL_END_NAMESPACE

#endif // CGAL_VALUE_TYPE_TRAITS_H

