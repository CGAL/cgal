// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : SoX
// File          : include/SoX/GAPS/macros.h
// SoX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file SoX/GAPS/macros.h
    \brief some basic macro declarations
*/

#ifndef SoX_GAPS_MACROS_H
#define SoX_GAPS_MACROS_H 1

#include <CGAL/config.h>

// TASK document 
#define SoX_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(Traits) \
 typedef typename Traits::Surface_3 Surface_3; \
 typedef typename Traits::Arrangement_traits_2 Arrangement_traits_2; \
 typedef typename Arrangement_traits_2::Curve_kernel_2 Curve_kernel_2; \
 typedef typename Traits::Z_at_xy_isolator Z_at_xy_isolator;\
 typedef typename Traits::Polynomial_3 Polynomial_3; \
 typedef typename Traits::Polynomial_2 Polynomial_2; \
 typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2; \
 typedef typename Traits::Point_2 Point_2; \
 typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2; \
// end #define SoX_SURFACE_Z_AT_XY_ISOLATOR_TRAITS_SNAP_TYPEDEFS(Traits)


#endif // SoX_GAPS_MACROS_H
// EOF
