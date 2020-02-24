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
// Author(s)     : Stefan Schirra


#ifndef CGAL_CH_EDDY_H
#define CGAL_CH_EDDY_H

#include <CGAL/license/Convex_hull_2.h>


#include <CGAL/basic.h>
#include <iterator>

namespace CGAL {

// same as |convex_hull_2(first,last,result)|. {\sc traits}: uses 
// |Traits::Point_2|, |Traits::Less_signed_distance_to_line_2|,
// |Traits::Left_turn_2|,  |Traits::Equal_2| and |Traits::Less_xy_2|.
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_eddy(InputIterator first, InputIterator last, 
        OutputIterator  result,
        const Traits& ch_traits);

template <class InputIterator, class OutputIterator>
inline
OutputIterator
ch_eddy(InputIterator first, InputIterator last, OutputIterator  result)
{
    typedef std::iterator_traits<InputIterator>   ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return ch_eddy( first, last, result, Kernel());
}

} //namespace CGAL

#include <CGAL/Convex_hull_2/ch_eddy_impl.h>

#endif // CGAL_CH_EDDY_H
