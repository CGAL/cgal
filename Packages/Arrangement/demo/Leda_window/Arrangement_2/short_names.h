// Copyright (c) 2003, 2004  Tel-Aviv University (Israel).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#include <CGAL/config.h>

#if defined(CGAL_CFG_LONGNAME_BUG)
//#define Quotient                        Qt
#define Homogeneous                     Hs
#define Cartesian                       Cn
#define Simple_cartesian                SC
#define Filtered_kernel                 FKl
#define Segment_2                       St
#define Point_2                         Pt2
#define Topological_map                 TM
#define Planar_map_2                    PMp
#define Arrangement_2                   Ar
#define I_HalfedgeDS_iterator           IPI
#define I_HalfedgeDS_const_iterator     IPCI
#define Arr_2_default_dcel              ADD
#define Arr_segment_traits_2            AST
#define Arr_circles_real_traits         ACRT
#define Arr_polyline_traits             APT
#define Arr_base_node                   ABN
#define Arr_2_halfedge_base             AHB
#define Arr_2_vertex_base               AVB
#define Arr_2_face_base                 AFB
#if ! defined(_MSC_VER)
#define bidirectional_iterator_tag      BIT
#endif
#define In_place_list_iterator          IPLI
#define In_place_list_const_iterator    IPLCI
#define allocator                       All
#define Pm_traits_wrap_2                PmTW
#define Td_X_trapezoid                  TXT
#define PL_X_curve_plus                 PXCP

#if defined(_MSC_VER)
// Has no effect, probably bug in MSVC
#pragma warning(disable:4503)
#endif

#endif
