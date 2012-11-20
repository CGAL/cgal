// Copyright (c) 1997  
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
// 
//
// Author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>

#ifndef CGAL_CIRCULATOR_BASES_H
#define CGAL_CIRCULATOR_BASES_H 1

#include <cstddef>
#include <iterator>

namespace CGAL {

struct Circulator_tag {};                   // any circulator.
struct Iterator_tag {};                     // any iterator.

// conversion operators instead of inheritance to avoid ambiguous
// bases. we have to repeat all possible conversion so we don't run
// into a multiple user-defined conversions, problem.

struct Forward_circulator_tag
  : public std::forward_iterator_tag 
{};

struct Bidirectional_circulator_tag
  : public std::bidirectional_iterator_tag 
{ 
  operator Forward_circulator_tag() const { return Forward_circulator_tag(); }
};

struct Random_access_circulator_tag
  : public std::random_access_iterator_tag
{ 
  operator Bidirectional_circulator_tag() const { return Bidirectional_circulator_tag(); }
  operator Forward_circulator_tag() const { return Forward_circulator_tag(); }
};

template <typename Tag, typename T, typename Distance = std::ptrdiff_t,
          /* size is so awkwardly placed to faciliate using the
           * default arguments from the derived classes */
          typename Size = std::size_t, typename Pointer = T*, 
          typename Reference = T&>
struct Circulator_base {
  typedef Tag       iterator_category;
  typedef T         value_type;
  typedef Distance  difference_type;
  typedef Pointer   pointer;
  typedef Reference reference;
  typedef Size      size_type;
};

template <class T, class Dist = std::ptrdiff_t, class Size = std::size_t>
struct Forward_circulator_base 
  : Circulator_base<Forward_circulator_tag, T, Dist, Size> {};

template <class T, class Dist = std::ptrdiff_t, class Size = std::size_t>
struct Bidirectional_circulator_base 
  : Circulator_base<Bidirectional_circulator_tag, T, Dist, Size> {};

template <class T, class Dist = std::ptrdiff_t, class Size = std::size_t>
struct Random_access_circulator_base
  : Circulator_base<Random_access_circulator_tag, T, Dist, Size> {};

// variant base classes
// ---------------------
template <class T, class Dist = std::ptrdiff_t, class Size = std::size_t>
class Forward_circulator_ptrbase         // forward circulator.
{
    protected:
        void* _ptr;
    public:
        typedef Forward_circulator_tag  iterator_category;
        typedef T                           value_type;
        typedef Dist                        difference_type;
        typedef Size                        size_type;
        typedef T*                          pointer;
        typedef T&                          reference;
        Forward_circulator_ptrbase()        : _ptr(NULL) {}
        Forward_circulator_ptrbase(void* p) : _ptr(p) {}
};
template <class T, class Dist = std::ptrdiff_t, class Size = std::size_t>
class Bidirectional_circulator_ptrbase   // bidirectional circulator.
{
    protected:
        void* _ptr;
    public:
        typedef Bidirectional_circulator_tag  iterator_category;
        typedef T                           value_type;
        typedef Dist                        difference_type;
        typedef Size                        size_type;
        typedef T*                          pointer;
        typedef T&                          reference;
        Bidirectional_circulator_ptrbase()        : _ptr(NULL) {}
        Bidirectional_circulator_ptrbase(void* p) : _ptr(p) {}
};
template <class T, class Dist = std::ptrdiff_t, class Size = std::size_t>
class Random_access_circulator_ptrbase   // random access circulator.
{
    protected:
        void* _ptr;
    public:
        typedef Random_access_circulator_tag iterator_category;
        typedef T                           value_type;
        typedef Dist                        difference_type;
        typedef Size                        size_type;
        typedef T*                          pointer;
        typedef T&                          reference;
        Random_access_circulator_ptrbase()        : _ptr(NULL) {}
        Random_access_circulator_ptrbase(void* p) : _ptr(p) {}
};

} //namespace CGAL

#endif // CGAL_CIRCULATOR_BASES_H //
// EOF //
