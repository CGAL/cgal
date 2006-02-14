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

#ifndef CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
#define CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Lazy_exact_nt.h>

#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpq.h>
#else
#  include <CGAL/Quotient.h>
#  include <CGAL/MP_Float.h>
#endif

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_USE_GMP
typedef Filtered_kernel<Simple_cartesian<Lazy_exact_nt<Gmpq > > >
        Exact_predicates_exact_constructions_kernel;
#else
typedef Filtered_kernel<Simple_cartesian<Lazy_exact_nt<Quotient<MP_Float> > > >
        Exact_predicates_exact_constructions_kernel;
#endif

CGAL_END_NAMESPACE

#endif // CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
