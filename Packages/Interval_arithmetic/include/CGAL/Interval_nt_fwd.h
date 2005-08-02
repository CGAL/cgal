// Copyright (c) 1998-2005  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_INTERVAL_NT_FWD_H
#define CGAL_INTERVAL_NT_FWD_H

// Forward declarations of functions based on Interval_nt, required by
// two-stage name lookup.

#include <CGAL/Uncertain.h>

CGAL_BEGIN_NAMESPACE

template <bool> class Interval_nt;

template <bool Protected>
double to_double (const Interval_nt<Protected> &);

template <bool Protected>
std::pair<double, double> to_interval (const Interval_nt<Protected> &);

template <bool Protected>
bool is_valid (const Interval_nt<Protected> &);

template <bool Protected>
bool is_finite (const Interval_nt<Protected> &);

template <bool Protected>
Interval_nt<Protected> sqrt (const Interval_nt<Protected> &);

template <bool Protected>
Interval_nt<Protected>
min (const Interval_nt<Protected> &, const Interval_nt<Protected> &);

template <bool Protected>
Interval_nt<Protected>
max (const Interval_nt<Protected> &, const Interval_nt<Protected> &);

template <bool Protected>
Interval_nt<Protected> square (const Interval_nt<Protected> &);

template <bool Protected>
Interval_nt<Protected> abs (const Interval_nt<Protected> &);

template <bool Protected>
Uncertain<Sign> sign (const Interval_nt<Protected> &);

template <bool Protected>
Uncertain<Comparison_result>
compare (const Interval_nt<Protected> &, const Interval_nt<Protected> &);

template <bool Protected>
Uncertain<bool> is_zero (const Interval_nt<Protected> & d);

template <bool Protected>
Uncertain<bool> is_one (const Interval_nt<Protected> & d);

template <bool Protected>
Uncertain<bool>
is_positive (const Interval_nt<Protected> & d);

template <bool Protected>
Uncertain<bool>
is_negative (const Interval_nt<Protected> & d);

CGAL_END_NAMESPACE

#endif // CGAL_INTERVAL_NT_FWD_H
