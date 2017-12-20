// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Ron Wein           <wein@post.tau.ac.il>
//                 Baruch Zukerman    <baruchzu@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_3_H
#define CGAL_ENVELOPE_3_H

#include <CGAL/license/Envelope_3.h>


#include <CGAL/Envelope_3/Envelope_diagram_on_surface_2.h>
#include <CGAL/Envelope_3/Envelope_divide_and_conquer_3.h>
#include <CGAL/Arr_accessor.h>

namespace CGAL {

/*!
 * Construct the lower envelope of a given set of surfaces.
 * \param begin An iterator for the first surface.
 * \param end A past-the-end iterator for the surfaces in the range.
 * \param min_diag Output: The minimization diagram.
 * \pre The value-type of InputIterator is Traits::Surface_3.
 */
template <typename InputIterator, typename GeomTraits,
          class TopTraits>
void lower_envelope_3 (InputIterator begin, InputIterator end,
                       Envelope_diagram_on_surface_2<GeomTraits, TopTraits>& 
                         min_diagram)
{
  typedef GeomTraits                                        Traits_3;
  typedef typename Envelope_diagram_on_surface_2<Traits_3, 
    TopTraits>::Arrangement                                 Envelope_diagram_2;
  typedef Envelope_divide_and_conquer_3<Traits_3, Envelope_diagram_2> 
                                                            Envelope_algorithm;
  Envelope_algorithm   env_alg (min_diagram.geometry_traits(), ENVELOPE_LOWER);
  env_alg.construct_lu_envelope (begin, end, min_diagram);
}

/*!
 * Construct the upper envelope of a given set of surfaces.
 * \param begin An iterator for the first surface.
 * \param end A past-the-end iterator for the surfaces in the range.
 * \param max_diag Output: The maximization diagram.
 * \pre The value-type of InputIterator is Traits::Surface_3.
 */
template <class InputIterator, class GeomTraits,
          class TopTraits>
void upper_envelope_3 (InputIterator begin, InputIterator end,
                       Envelope_diagram_on_surface_2<GeomTraits, TopTraits>& 
                         max_diagram)
{
  typedef GeomTraits                                         Traits_3;
  typedef typename Envelope_diagram_on_surface_2<
  Traits_3, TopTraits>::Arrangement                          Envelope_diagram_2;
  typedef Envelope_divide_and_conquer_3<Traits_3, Envelope_diagram_2> 
                                                            Envelope_algorithm;

  Envelope_algorithm   env_alg (max_diagram.geometry_traits(), ENVELOPE_UPPER);
  env_alg.construct_lu_envelope (begin, end, max_diagram);
}


/*!
 * Construct the lower envelope of a given set of xy_monotone surfaces.
 * \param begin An iterator for the first xy-monotone surface.
 * \param end A past-the-end iterator for the xy_monotone surfaces in the range.
 * \param min_diag Output: The minimization diagram.
 * \pre The value-type of InputIterator is GeomTraits::Surface_3.
 */
template <class InputIterator, class GeomTraits,
          class TopTraits>
void
lower_envelope_xy_monotone_3 (InputIterator begin, InputIterator end,
                              Envelope_diagram_on_surface_2<GeomTraits, 
                              TopTraits>& min_diagram)
{
  typedef GeomTraits                                         Traits_3;
  typedef typename Envelope_diagram_on_surface_2<
  Traits_3, TopTraits>::Arrangement                          Envelope_diagram_2;
  typedef Envelope_divide_and_conquer_3<Traits_3, Envelope_diagram_2> 
                                                            Envelope_algorithm;
  Envelope_algorithm   env_alg (min_diagram.geometry_traits(), ENVELOPE_LOWER);
  env_alg.construct_envelope_xy_monotone (begin, end, min_diagram);
}


/*!
 * Construct the upper envelope of a given set of xy_monotone surfaces.
 * \param begin An iterator for the first xy_monotone surface.
 * \param end A past-the-end iterator for the xy_monotone surfaces in the range.
 * \param max_diag Output: The maximization diagram.
 * \pre The value-type of InputIterator is Traits::Surface_3.
 */
template <class InputIterator, class GeomTraits,
          class TopTraits>
void
upper_envelope_xy_monotone_3 (InputIterator begin, InputIterator end,
                              Envelope_diagram_on_surface_2<GeomTraits,
                              TopTraits>& max_diagram)
{
  typedef GeomTraits                                         Traits_3;
  typedef typename Envelope_diagram_on_surface_2<
  Traits_3, TopTraits>::Arrangement                          Envelope_diagram_2;
  typedef Envelope_divide_and_conquer_3<Traits_3, Envelope_diagram_2> 
                                                            Envelope_algorithm;

  Envelope_algorithm   env_alg (max_diagram.geometry_traits(), ENVELOPE_UPPER);
  env_alg.construct_envelope_xy_monotone (begin, end, max_diagram);

  return;
}


} //namespace CGAL

#endif
