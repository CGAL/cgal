// ============================================================================
//
// Copyright (c) 1997,1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-wip $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/workaround_casts.h
// chapter       : $CGAL_Chapter: Configuration $
// package       : $CGAL_Package: Workarounds $
//
// source        : web/workarounds.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sven Schönherr <sven@inf.fu-berlin.de>
//
// coordinator   : Utrecht University (Wieger Wesselink <wieger@cs.ruu.nl>)
//
// implementation: Workarounds for C++-style casts
// ============================================================================

#ifndef CGAL_WORKAROUND_CASTS_H
#define CGAL_WORKAROUND_CASTS_H 1


// workaround for C++-style casts

#if defined( CGAL_CFG_NO_DYNAMIC_CAST)
#  define  CGAL_dynamic_cast(type,expr)  (type)( expr)
#else
#  define  CGAL_dynamic_cast(type,expr)  dynamic_cast< type >(expr)
#endif // CGAL_CFG_NO_DYNAMIC_CAST

#endif // CGAL_WORKAROUND_CASTS_H

// ===== EOF ==================================================================
