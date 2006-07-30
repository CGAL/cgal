// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source: $
// $Revision$ $Date$
// $Name:  $
//
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_2_H
#define CGAL_ENVELOPE_2_H

/*! \file
 * Global functions fro computing lower and upper envelopes of curves in the
 * plane.
 */

#include <Envelope_2/Env_divide_and_conquer_2.h>

CGAL_BEGIN_NAMESPACE

/*!
 * Compute the lower envelope of a range of curves.
 * \param traits A traits class instance.
 * \param begin An iterator for the first curve.
 * \param end A past-the-end iterator for the curves.
 * \param diag Output: The minimization diagram.
 * \pre The value-type of the iterator is Traits::Curve_2.
 */
template <class Traits, class InputIterator, class EnvelopeDiagram>
void lower_envelope_2 (Traits& traits,
                       InputIterator begin, InputIterator end,
                       EnvelopeDiagram& diag)
{
  typedef Envelope_divide_and_conquer_2<Traits, Diagram>   Envelope_2;

  Envelope_2      env (&traits);

  env.insert_curves (begin, end,
                     typename Envelope_2::LOWER,
                     diag);

  return;
}

/*!
 * Compute the upper envelope of a range of curves.
 * \param traits A traits class instance.
 * \param begin An iterator for the first curve.
 * \param end A past-the-end iterator for the curves.
 * \param diag Output: The maximization diagram.
 * \pre The value-type of the iterator is Traits::Curve_2.
 */
template <class Traits, class InputIterator, class EnvelopeDiagram>
void upper_envelope_2 (Traits& traits,
                       InputIterator begin, InputIterator end,
                       EnvelopeDiagram& diag)
{
  typedef Envelope_divide_and_conquer_2<Traits, Diagram>   Envelope_2;

  Envelope_2      env (&traits);

  env.insert_curves (begin, end,
                     typename Envelope_2::UPPER,
                     diag);

  return;
}

/*!
 * Compute the lower envelope of a range of x-monotone curves.
 * \param traits A traits class instance.
 * \param begin An iterator for the first x-monotone curve.
 * \param end A past-the-end iterator for the x-monotone curves.
 * \param diag Output: The minimization diagram.
 * \pre The value-type of the iterator is Traits::X_monotone_curve_2.
 */
template <class Traits, class InputIterator, class EnvelopeDiagram>
void lower_envelope_x_monotone_2 (Traits& traits,
                                  InputIterator begin, InputIterator end,
                                  EnvelopeDiagram& diag)
{
  typedef Envelope_divide_and_conquer_2<Traits, Diagram>   Envelope_2;

  Envelope_2      env (&traits);

  env.insert_x_monotone_curves (begin, end,
                                typename Envelope_2::LOWER,
                                diag);

  return;
}

/*!
 * Compute the upper envelope of a range of x-monotone curves.
 * \param traits A traits class instance.
 * \param begin An iterator for the first x-monotone curve.
 * \param end A past-the-end iterator for the x-monotone curves.
 * \param diag Output: The maximization diagram.
 * \pre The value-type of the iterator is Traits::X_monotone_curve_2.
 */
template <class Traits, class InputIterator, class EnvelopeDiagram>
void upper_envelope_x_monotone_2 (Traits& traits,
                                  InputIterator begin, InputIterator end,
                                  EnvelopeDiagram& diag)
{
  typedef Envelope_divide_and_conquer_2<Traits, Diagram>   Envelope_2;

  Envelope_2      env (&traits);

  env.insert_x_monotone_curves (begin, end,
                                typename Envelope_2::UPPER,
                                diag);

  return;
}

CGAL_END_NAMESPACE

#endif
