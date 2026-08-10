// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ron Wein   <wein@post.tau.ac.il>
//            Efi Fogel  <efifogel@gmail.com>

#ifndef CGAL_ENVELOPE_2_H
#define CGAL_ENVELOPE_2_H

#include <CGAL/license/Envelope_2.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/envelope_2.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/Envelope_2/envelope_2>"
#include <CGAL/Installation/internal/deprecation_warning.h>

#if (defined __GNUC__)
  #if !(defined __STRICT_ANSI__)
  #warning "envelope_2.h is DEPRECATED, please include Envelope_2/envelope_2.h instead."
  #endif
#elif (defined _MSC_VER)
  #pragma message("envelope_2.h is DEPRECATED, please include Envelope_2/envelope_2.h instead")
#endif

#include <CGAL/Envelope_2/envelope_2.h>

/*! \file
 * Deprecated global functions for computing lower and upper envelopes of curves in the plane.
 */

namespace CGAL {

/*! computes the lower envelope of a range of curves.
 * \param begin An iterator for the first curve.
 * \param end A past-the-end iterator for the curves.
 * \param diag Output: The minimization diagram.
 * \pre The value-type of the iterator is Traits::Curve_2.
 */
template <typename InputIterator, typename EnvelopeDiagram>
CGAL_DEPRECATED void lower_envelope_2(InputIterator begin, InputIterator end, EnvelopeDiagram& diag)
{ Envelope_2::lower_envelope_2(begin, end, diag); }

/*! computes the upper envelope of a range of curves.
 * \param begin An iterator for the first curve.
 * \param end A past-the-end iterator for the curves.
 * \param diag Output: The maximization diagram.
 * \pre The value-type of the iterator is Traits::Curve_2.
 */
template <typename InputIterator, typename EnvelopeDiagram>
CGAL_DEPRECATED void upper_envelope_2(InputIterator begin, InputIterator end, EnvelopeDiagram& diag)
{ Envelope_2::upper_envelope_2(begin, end, diag); }

/*! computes the lower envelope of a range of x-monotone curves.
 * \param begin An iterator for the first x-monotone curve.
 * \param end A past-the-end iterator for the x-monotone curves.
 * \param diag Output: The minimization diagram.
 * \pre The value-type of the iterator is Traits::X_monotone_curve_2.
 */
template <typename InputIterator, typename EnvelopeDiagram>
CGAL_DEPRECATED void lower_envelope_x_monotone_2(InputIterator begin, InputIterator end, EnvelopeDiagram& diag)
{ Envelope_2::lower_envelope_x_monotone_2(begin, end, diag); }

/*! computes the lower envelope of a range of x-monotone curves.
 * Compute the lower envelope of a range of x-monotone curves.
 * \param begin An iterator for the first x-monotone curve.
 * \param end A past-the-end iterator for the x-monotone curves.
 * \param diag Output: The minimization diagram.
 * \param traits The arrangement traits responsible for the x-monotone curves.
 * \pre The value-type of the iterator is Traits::X_monotone_curve_2.
 */
template <typename InputIterator, typename EnvelopeDiagram, typename Traits>
CGAL_DEPRECATED void lower_envelope_x_monotone_2(InputIterator begin, InputIterator end, EnvelopeDiagram& diag,
                                                 const Traits& traits)
{ Envelope_2::lower_envelope_x_monotone_2(begin, end, diag, traits); }

/*! computes the upper envelope of a range of x-monotone curves.
 * \param begin An iterator for the first x-monotone curve.
 * \param end A past-the-end iterator for the x-monotone curves.
 * \param diag Output: The maximization diagram.
 * \pre The value-type of the iterator is Traits::X_monotone_curve_2.
 */
template <typename InputIterator, typename EnvelopeDiagram>
CGAL_DEPRECATED void upper_envelope_x_monotone_2(InputIterator begin, InputIterator end, EnvelopeDiagram& diag)
{ Envelope_2::upper_envelope_x_monotone_2(begin, end, diag); }

/*! computes the upper envelope of a range of x-monotone curves.
 * \param begin An iterator for the first x-monotone curve.
 * \param end A past-the-end iterator for the x-monotone curves.
 * \param diag Output: The maximization diagram.
 * \param traits The arrangement traits responsible for the x-monotone curves.
 * \pre The value-type of the iterator is Traits::X_monotone_curve_2.
 */
template <typename InputIterator, typename EnvelopeDiagram, typename Traits>
CGAL_DEPRECATED void upper_envelope_x_monotone_2(InputIterator begin, InputIterator end, EnvelopeDiagram& diag,
                                                 const Traits& traits)
{ Envelope_2::upper_envelope_x_monotone_2(begin, end, diag, traits); }

} // namespace CGAL

#endif
