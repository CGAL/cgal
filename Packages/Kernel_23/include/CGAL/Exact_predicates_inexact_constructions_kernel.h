// Copyright (c) 2003  Utrecht University (The Netherlands),
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
// Author(s)     : Menelaos Karavelas, Sylvain Pion

#ifndef CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
#define CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H

#include <CGAL/Simple_cartesian.h>
// We don't use Filtered_kernel at the moment because it's slower
// that Filtered_exact.
// #include <CGAL/Filtered_kernel.h>

#ifdef CGAL_CFG_MATCHING_BUG_2
#  define CGAL_IA_CT        double
#  define CGAL_IA_PROTECTED true
#  define CGAL_IA_CACHE     No_Filter_Cache
#  define CGAL_IA_ET        CGAL::MP_Float
#endif

#include <CGAL/MP_Float.h>
#include <CGAL/Filtered_exact.h>

CGAL_BEGIN_NAMESPACE

#if 0
typedef Filtered_kernel< Simple_cartesian<double> >
        Exact_predicates_inexact_constructions_kernel;
#else
typedef Simple_cartesian<Filtered_exact<double, MP_Float> >
        Exact_predicates_inexact_constructions_kernel;
#endif

CGAL_END_NAMESPACE

#endif // CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
