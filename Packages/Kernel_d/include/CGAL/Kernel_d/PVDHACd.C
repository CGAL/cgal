// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Kernel_d/PVDHACd.C
// package       : Kernel_d
// chapter       : Basic
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Susan Hert <hert@mpi-sb.mpg.de>
//
// implementation: implementation inclusion
// ============================================================================

#if defined(CGAL_POINTCD_H) && defined(CGAL_VECTORCD_H) && \
    defined(CGAL_DIRECTIONCD_H) && defined(CGAL_HYPERPLANECD_H) && \
    defined(CGAL_AFF_TRANSFORMATIONCD_H) && !defined(CGAL_PVDHACD_C)
#define CGAL_PVDHACD_C

#include <CGAL/Kernel_d/PointCd.C>
#include <CGAL/Kernel_d/VectorCd.C>
#include <CGAL/Kernel_d/DirectionCd.C> 
#include <CGAL/Kernel_d/HyperplaneCd.C>
#include <CGAL/Kernel_d/Aff_transformationCd.C>

#endif //CGAL_PVDHACD_C


