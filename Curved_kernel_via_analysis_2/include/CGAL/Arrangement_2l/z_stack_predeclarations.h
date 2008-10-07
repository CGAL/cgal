// Copyright (c) 2007-2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>

#ifndef CGAL_ARRANGEMENT_2l_Z_STACK_PREDECLARATIONS_H
#define CGAL_ARRANGEMENT_2l_Z_STACK_PREDECLARATIONS_H 1

/*!\file include/CGAL/Arrangement_2l/z_stack_predeclarations.h
 * \brief contains all required pre-declarations
 */

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

#endif // CGAL_ARRANGEMENT_2l_Z_STACK_PREDECLARATIONS_H
// EOF
