// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>

#ifndef CGAL_ARITHMETIC_FILTER_DISPATCH_H
#define CGAL_ARITHMETIC_FILTER_DISPATCH_H

#if !defined( CGAL_ARITHMETIC_FILTER_PREDICATES_BUILTIN_H )
#include <CGAL/Arithmetic_filter/predicates/builtin.h>
#endif

#if defined( CGAL_PREDICATES_SIGN_OF_DETERMINANT_H ) && \
       !defined( CGAL_ARITHMETIC_FILTER_PREDICATES_SIGN_OF_DETERMINANT_H )
#include <CGAL/Arithmetic_filter/predicates/sign_of_determinant.h>
#endif

#if defined( CGAL_PREDICATES_KERNEL_FTC2_H ) && \
   !defined( CGAL_ARITHMETIC_FILTER_PREDICATES_KERNEL_FTC2_H )
#include <CGAL/Arithmetic_filter/predicates/kernel_ftC2.h>
#endif

#if defined( CGAL_PREDICATES_KERNEL_FTC3_H ) && \
       !defined( CGAL_ARITHMETIC_FILTER_PREDICATES_KERNEL_FTC3_H )
#include <CGAL/Arithmetic_filter/predicates/kernel_ftC3.h>
#endif

#if defined( CGAL_REGULAR_TRIANGULATION_FTC2_H ) && \
       !defined( CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_FTC2_H )
#include <CGAL/Arithmetic_filter/predicates/Regular_triangulation_ftC2.h>
#endif

    /*
#if defined( CGAL_REGULAR_TRIANGULATION_RTH2_H ) && \
       !defined( CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_RTH2_H )
#include <CGAL/Arithmetic_filter/predicates/Regular_triangulation_rtH2.h>
#endif
*/

#if defined( CGAL_REGULAR_TRIANGULATION_FTC3_H ) && \
       !defined( CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_FTC3_H )
#include <CGAL/Arithmetic_filter/predicates/Regular_triangulation_ftC3.h>
#endif

    /*
#if defined( CGAL_REGULAR_TRIANGULATION_RTH3_H ) && \
       !defined( CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_RTH3_H )
#include <CGAL/Arithmetic_filter/predicates/Regular_triangulation_rtH3.h>
#endif
*/

#endif // CGAL_ARITHMETIC_FILTER_DISPATCH_H
