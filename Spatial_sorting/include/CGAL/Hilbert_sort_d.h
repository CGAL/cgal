// Copyright (c) 2011  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Olivier Devillers

#ifndef CGAL_HILBERT_SORT_d_H
#define CGAL_HILBERT_SORT_d_H

#include <CGAL/Hilbert_policy_tags.h>
#include <CGAL/Hilbert_sort_median_d.h>
#include <CGAL/Hilbert_sort_middle_d.h>

namespace CGAL {

template <class K,  class Hilbert_policy >
class Hilbert_sort_d;

template <class K>
class Hilbert_sort_d<K, Hilbert_sort_median_policy >
    : public Hilbert_sort_median_d<K>
{
public:
  Hilbert_sort_d (const K &k=K() , std::ptrdiff_t limit=1 )
    : Hilbert_sort_median_d<K> (k,limit)
  {}
};

template <class K>
class Hilbert_sort_d<K, Hilbert_sort_middle_policy >
    : public Hilbert_sort_middle_d<K>
{
public:
  Hilbert_sort_d (const K &k=K() , std::ptrdiff_t limit=1 )
    : Hilbert_sort_middle_d<K> (k,limit)
  {}
};

} // namespace CGAL

#endif//CGAL_HILBERT_SORT_d_H
