// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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

struct Forward_circulator_tag
    : public std::forward_iterator_tag {};
struct Bidirectional_circulator_tag
    : public std::bidirectional_iterator_tag {};
struct Random_access_circulator_tag
    : public std::random_access_iterator_tag {};
template <class T, class Dist = std::ptrdiff_t, class Size = std::size_t>
struct Forward_circulator_base {
    typedef T                            value_type;
    typedef Dist                         difference_type;
    typedef Size                         size_type;
    typedef T*                           pointer;
    typedef T&                           reference;
    typedef Forward_circulator_tag       iterator_category;
};
template <class T, class Dist = std::ptrdiff_t, class Size = std::size_t>
struct Bidirectional_circulator_base {
    typedef T                            value_type;
    typedef Dist                         difference_type;
    typedef Size                         size_type;
    typedef T*                           pointer;
    typedef T&                           reference;
    typedef Bidirectional_circulator_tag iterator_category;
};
template <class T, class Dist = std::ptrdiff_t, class Size = std::size_t>
struct Random_access_circulator_base {
    typedef T                            value_type;
    typedef Dist                         difference_type;
    typedef Size                         size_type;
    typedef T*                           pointer;
    typedef T&                           reference;
    typedef Random_access_circulator_tag iterator_category;
};
template < class Category,
           class T,
           class Distance  = std::ptrdiff_t,
           class Size      = std::size_t,
           class Pointer   = T*,
           class Reference = T&>
struct Circulator_base {
    typedef Category  iterator_category;
    typedef T         value_type;
    typedef Distance  difference_type;
    typedef Size      size_type;
    typedef Pointer   pointer;
    typedef Reference reference;
};

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
