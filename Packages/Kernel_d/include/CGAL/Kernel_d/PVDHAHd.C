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
// file          : include/CGAL/Kernel_d/PVDHAHd.C
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

#if defined(CGAL_POINTHD_H) && defined(CGAL_VECTORHD_H) && \
    defined(CGAL_DIRECTIONHD_H) && defined(CGAL_HYPERPLANEHD_H) && \
    defined(CGAL_AFF_TRANSFORMATIONHD_H) && !defined(CGAL_PVDHAHD_C)
#define CGAL_PVDHAHD_C

#include <CGAL/Kernel_d/PointHd.C>
#include <CGAL/Kernel_d/VectorHd.C>
#include <CGAL/Kernel_d/DirectionHd.C> 
#include <CGAL/Kernel_d/HyperplaneHd.C>
#include <CGAL/Kernel_d/Aff_transformationHd.C>

#endif //CGAL_PVDHAHD_C


