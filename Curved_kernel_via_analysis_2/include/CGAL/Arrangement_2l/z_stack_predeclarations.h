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
// File          : include/SoX/GAPS/z_stack_predclarations.h
// SoX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file SoX/GAPS/z_stack_predeclarations.h
    \brief contains all required pre-declarations
*/

#ifndef SoX_GAPS_Z_STACK_PREDECLARATIONS_H
#define SoX_GAPS_Z_STACK_PREDECLARATIONS_H 1

#include <CGAL/config.h>

CGAL_BEGIN_NAMESPACE

// pre-declaration
template < class SurfaceZAtXyIsolatorTraits >
class Restricted_cad_3;

// pre-declaration
template < class SurfaceZAtXyIsolatorTraits >
class P_dcel_info;

// pre-declaration
template < class SurfaceZAtXyIsolatorTraits, class DcelData >
class Z_stack;

// pre-declaration
template < class SurfaceZAtXyIsolatorTraits >
class Create_restricted_cad_3;

// pre-declaration
template < class SurfaceZAtXyIsolatorTraits >
class Overlay_restricted_cad_3;

// pre-declaration
template < class Arrangement_2_ >
class Arr_p_dcel_info_overlay_traits;

// pre-declaration
template < class SurfaceZAtXyIsolatorTraits, class Rep_ >
class Surface_pair_3;

namespace CGALi {

// pre-declaration
template < class SurfaceZAtXyIsolatorTraits >
class Restricted_cad_3_cache;

// pre-declaration
template < class SurfaceZAtXyIsolatorTraits, class DcelData >
class Z_cell;

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // SoX_GAPS_Z_STACK_PREDECLARATIONS_H
// EOF
