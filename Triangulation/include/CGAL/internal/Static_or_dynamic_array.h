// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)    : Samuel Hornus

#ifndef CGAL_INTERNAL_STATIC_OR_DYNAMIC_ARRAY_H
#define CGAL_INTERNAL_STATIC_OR_DYNAMIC_ARRAY_H

#include <CGAL/license/Triangulation.h>


#include <CGAL/Compact_container.h>
#include <CGAL/Dimension.h>
#include <CGAL/array.h>
#include <vector>

namespace CGAL {

namespace internal {

// Utility for adding one to an Dimension_tag:

template<typename D>
struct Dimen_plus_one;

template<>
struct Dimen_plus_one<Dynamic_dimension_tag>
{
    typedef Dynamic_dimension_tag type;
};

template<int D>
struct Dimen_plus_one<Dimension_tag<D> >
{
    typedef Dimension_tag<D+1> type;
};

// A SMALL CONTAINER UTILITY FOR DYNAMIC/STATIC MEMORY MANAGEMENT

// stores an array of static or dynamic size, depending on template parameter <B>.

template< typename Containee, typename D, bool WithCompactContainerHelper = false>
    struct S_or_D_array; // S = static, D = dynamic

// The case of static size:
template< typename Containee, int D, bool WithCompactContainerHelper >
struct S_or_D_array< Containee, Dimension_tag< D >, WithCompactContainerHelper >
: public std::array<Containee, D>
{
    typedef std::array<Containee, D> Base;
    S_or_D_array(const int)
    : Base()
    {}
    S_or_D_array(const int, const Containee & c)
    : Base()
    {
        assign(c);
    }
    void*   for_compact_container() const
    {
        return (*this)[0].for_compact_container();
    }
    void    for_compact_container(void *p)
    {
        (*this)[0].for_compact_container(p);
    }
};

// The case of dynamic size
template< typename Containee >
struct S_or_D_array< Containee, Dynamic_dimension_tag, false >
: public std::vector<Containee>
{
    typedef std::vector<Containee> Base;
    // TODO: maybe we should use some "small-vector-optimized" class.
    S_or_D_array(const int d)
    : Base(d)
    {}
    S_or_D_array(const int d, const Containee & c)
    : Base(d, c)
    {}
};

// The case of dynamic size with for_compact_container
template< typename Containee >
struct S_or_D_array< Containee, Dynamic_dimension_tag, true >
: public std::vector<Containee>
{
    typedef std::vector<Containee> Base;
    S_or_D_array(const int d)
    : Base(d), fcc_(nullptr)
    {}
    S_or_D_array(const int d, const Containee & c)
    : Base(d, c), fcc_(nullptr)
    {}
    void* fcc_;
    void*   for_compact_container() const { return fcc_; }
    void    for_compact_container(void* p)       { fcc_ = p; }
};

} // end of namespace internal

} // end of namespace CGAL

#endif // CGAL_INTERNAL_STATIC_OR_DYNAMIC_ARRAY_H
