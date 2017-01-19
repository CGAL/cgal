// Copyright (c) 1999  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
