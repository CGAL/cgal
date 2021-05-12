// Copyright (c) 1999  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner


#ifndef CGAL_CH_MELKMAN_H
#define CGAL_CH_MELKMAN_H

#include <CGAL/license/Convex_hull_2.h>


#include <CGAL/basic.h>
#include <iterator>

namespace CGAL {

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_melkman( InputIterator first, InputIterator last,
            OutputIterator result, const Traits& ch_traits);


template <class InputIterator, class OutputIterator>
OutputIterator
ch_melkman( InputIterator first, InputIterator last,  OutputIterator result)
{
    typedef std::iterator_traits<InputIterator>   ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return ch_melkman( first, last, result, Kernel());
}

} //namespace CGAL

#include <CGAL/Convex_hull_2/ch_melkman_impl.h>

#endif // CGAL_CH_MELKMAN_H
