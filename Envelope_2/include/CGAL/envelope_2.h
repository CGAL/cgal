// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_2_H
#define CGAL_ENVELOPE_2_H

#include <CGAL/license/Envelope_2.h>


/*! \file
 * Global functions for computing lower and upper envelopes of curves in the
 * plane.
 */

#include <CGAL/Envelope_2/Env_divide_and_conquer_2.h>

namespace CGAL {

/*!
 * Compute the lower envelope of a range of curves.
 * \param begin An iterator for the first curve.
 * \param end A past-the-end iterator for the curves.
 * \param diag Output: The minimization diagram.
 * \pre The value-type of the iterator is Traits::Curve_2.
 */
template <class InputIterator, class EnvelopeDiagram>
void lower_envelope_2 (InputIterator begin, InputIterator end,
                       EnvelopeDiagram& diag)
{
  typedef typename EnvelopeDiagram::Traits_2                       Traits_2;
  typedef Envelope_divide_and_conquer_2<Traits_2, EnvelopeDiagram> Envelope_2;

  Envelope_2      env;

  env.insert_curves (begin, end,
                     true,         // Lower envelope.
                     diag);

  return;
}

/*!
 * Compute the upper envelope of a range of curves.
 * \param begin An iterator for the first curve.
 * \param end A past-the-end iterator for the curves.
 * \param diag Output: The maximization diagram.
 * \pre The value-type of the iterator is Traits::Curve_2.
 */
template <class InputIterator, class EnvelopeDiagram>
void upper_envelope_2 (InputIterator begin, InputIterator end,
                       EnvelopeDiagram& diag)
{
  typedef typename EnvelopeDiagram::Traits_2                       Traits_2;
  typedef Envelope_divide_and_conquer_2<Traits_2, EnvelopeDiagram> Envelope_2;

  Envelope_2      env;

  env.insert_curves (begin, end,
                     false,         // Upper envelope.
                     diag);

  return;
}

/*!
 * Compute the lower envelope of a range of x-monotone curves.
 * \param begin An iterator for the first x-monotone curve.
 * \param end A past-the-end iterator for the x-monotone curves.
 * \param diag Output: The minimization diagram.
 * \pre The value-type of the iterator is Traits::X_monotone_curve_2.
 */
template <class InputIterator, class EnvelopeDiagram>
void lower_envelope_x_monotone_2 (InputIterator begin, InputIterator end,
                                  EnvelopeDiagram& diag)
{
  typedef typename EnvelopeDiagram::Traits_2                       Traits_2;
  typedef Envelope_divide_and_conquer_2<Traits_2, EnvelopeDiagram> Envelope_2;

  Envelope_2      env;

  env.insert_x_monotone_curves (begin, end,
                                true,              // Lower envelope.
                                diag);

  return;
}

/*!
 * Compute the lower envelope of a range of x-monotone curves.
 * \param begin An iterator for the first x-monotone curve.
 * \param end A past-the-end iterator for the x-monotone curves.
 * \param diag Output: The minimization diagram.
 * \param traits The arrangement traits responsible for the x-monotone curves.
 * \pre The value-type of the iterator is Traits::X_monotone_curve_2.
 */
template <class InputIterator, class EnvelopeDiagram, class Traits>
void lower_envelope_x_monotone_2 (InputIterator begin, InputIterator end,
                                  EnvelopeDiagram& diag, const Traits& traits)
{
  typedef typename EnvelopeDiagram::Traits_2                       Traits_2;
  typedef Envelope_divide_and_conquer_2<Traits_2, EnvelopeDiagram> Envelope_2;

  Envelope_2      env{&traits};

  env.insert_x_monotone_curves (begin, end,
                                true,              // Lower envelope.
                                diag);

  return;
}

/*!
 * Compute the upper envelope of a range of x-monotone curves.
 * \param begin An iterator for the first x-monotone curve.
 * \param end A past-the-end iterator for the x-monotone curves.
 * \param diag Output: The maximization diagram.
 * \pre The value-type of the iterator is Traits::X_monotone_curve_2.
 */
template <class InputIterator, class EnvelopeDiagram>
void upper_envelope_x_monotone_2 (InputIterator begin, InputIterator end,
                                  EnvelopeDiagram& diag)
{
  typedef typename EnvelopeDiagram::Traits_2                       Traits_2;
  typedef Envelope_divide_and_conquer_2<Traits_2, EnvelopeDiagram> Envelope_2;

  Envelope_2      env;

  env.insert_x_monotone_curves (begin, end,
                                false,          // Upper envelope.
                                diag);

  return;
}

/*!
 * Compute the upper envelope of a range of x-monotone curves.
 * \param begin An iterator for the first x-monotone curve.
 * \param end A past-the-end iterator for the x-monotone curves.
 * \param diag Output: The maximization diagram.
 * \param traits The arrangement traits responsible for the x-monotone curves.
 * \pre The value-type of the iterator is Traits::X_monotone_curve_2.
 */
template <class InputIterator, class EnvelopeDiagram, class Traits>
void upper_envelope_x_monotone_2 (InputIterator begin, InputIterator end,
                                  EnvelopeDiagram& diag, const Traits& traits)
{
  typedef typename EnvelopeDiagram::Traits_2                       Traits_2;
  typedef Envelope_divide_and_conquer_2<Traits_2, EnvelopeDiagram> Envelope_2;

  Envelope_2      env{&traits};

  env.insert_x_monotone_curves (begin, end,
                                false,          // Upper envelope.
                                diag);

  return;
}

} //namespace CGAL

#endif
