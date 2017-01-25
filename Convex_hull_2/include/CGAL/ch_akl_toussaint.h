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


#ifndef CGAL_CH_AKL_TOUSSAINT_H
#define CGAL_CH_AKL_TOUSSAINT_H

#include <CGAL/license/Convex_hull_2.h>


#include <CGAL/basic.h>
#include <iterator>

namespace CGAL {

// same as |convex_hull_2(first,last,result)|.
// {\sc traits}: operates on |Traits::Point_2| using |Traits::Less_xy_2|, 
// |Traits::Less_yx_2|, |Traits::Equal_2| and |Traits::Left_turn_2|.
template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
ch_akl_toussaint(ForwardIterator first, ForwardIterator last, 
                 OutputIterator  result,
                 const Traits&   ch_traits);

template <class ForwardIterator, class OutputIterator>
inline
OutputIterator
ch_akl_toussaint(ForwardIterator first, ForwardIterator last, 
                 OutputIterator  result)
{
    typedef std::iterator_traits<ForwardIterator> ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return ch_akl_toussaint( first, last, result, Kernel()); 
}

} //namespace CGAL

#include <CGAL/Convex_hull_2/ch_akl_toussaint_impl.h>

#endif // CGAL_CH_AKL_TOUSSAINT_H
