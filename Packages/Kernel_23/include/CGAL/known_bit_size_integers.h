// Copyright (c) 1999,2001  Utrecht University (The Netherlands),
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

#ifndef CGAL_KNOWN_BIT_SIZE_INTEGERS_H
#define CGAL_KNOWN_BIT_SIZE_INTEGERS_H

CGAL_BEGIN_NAMESPACE

#if (defined(__sparc__) || defined(__sparc) || defined(sparc)) || \
    (defined(__sgi__)   || defined(__sgi)   || defined(sgi)) || \
    (defined(__i386__)  || defined(__i386)  || defined(i386)) || \
    (defined(__ia64__)  || defined(__ia64)  || defined(ia64)) || \
    (defined(__alpha__) || defined(__alpha) || defined(alpha)) || \
    (defined(__ppc__)   || defined(__ppc)   || defined(ppc)) || \
    (defined(__powerpc__) || defined(__powerpc) || defined(powerpc))
    typedef  signed char             Integer8;
    typedef  short                   Integer16;
    typedef  int                     Integer32;
    typedef  unsigned char           UInteger8;
    typedef  unsigned short          UInteger16;
    typedef  unsigned int            UInteger32;
    // See long_long.h for Integer64.
#  ifdef __ia64__
    typedef long                     Integer64;
    typedef unsigned long            UInteger64;
#    define CGAL_HAS_INTEGER64
#  endif
#else
#  if defined(__BORLANDC__)
    typedef  __int8                  Integer8;
    typedef  __int16                 Integer16;
    typedef  __int32                 Integer32;
    typedef  __int64                 Integer64;
    typedef  unsigned __int8         UInteger8;
    typedef  unsigned __int16        UInteger16;
    typedef  unsigned __int32        UInteger32;
    typedef  unsigned __int64        UInteger64;
#define CGAL_HAS_INTEGER64
#  else
#  if defined(_MSC_VER)
    typedef  signed char             Integer8;
    typedef  short                   Integer16;
    typedef  int                     Integer32;
    typedef  __int64                 Integer64;
    typedef  unsigned char           UInteger8;
    typedef  unsigned short          UInteger16;
    typedef  unsigned int            UInteger32;
    typedef  unsigned __int64        UInteger64;
#define CGAL_HAS_INTEGER64
#  else
#    error "patch this"
#  endif
#  endif
#endif

CGAL_END_NAMESPACE

#endif // CGAL_KNOWN_BIT_SIZE_INTEGERS_H
