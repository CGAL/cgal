// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// Author(s)     : Lutz Kettner
//                 Stefan Schirra
 

#ifndef CGAL_BASIC_H
#define CGAL_BASIC_H

#include <CGAL/config.h>

#include <iostream>
#include <cstdlib>


// Big endian or little endian machine.
// ====================================
#ifdef CGAL_CFG_NO_BIG_ENDIAN
#  define CGAL_LITTLE_ENDIAN 1
#else
#  define CGAL_BIG_ENDIAN 1
#endif

#ifndef CGAL_USE_LEDA
#  define CGAL_USE_CGAL_WINDOW
#endif

// Symbolic constants to tailor inlining. Inlining Policy.
// =======================================================
#ifndef CGAL_MEDIUM_INLINE
#  define CGAL_MEDIUM_INLINE inline
#endif

#ifndef CGAL_LARGE_INLINE
#  define CGAL_LARGE_INLINE
#endif

#ifndef CGAL_HUGE_INLINE
#  define CGAL_HUGE_INLINE
#endif

#include <CGAL/assertions.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/Object.h>
#include <CGAL/enum.h>
#include <CGAL/tags.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/IO/io.h>
#include <CGAL/Handle.h> // This should be removed ASAP.
#include <CGAL/kernel_basic.h>
#include <CGAL/known_bit_size_integers.h>

#endif // CGAL_BASIC_H
