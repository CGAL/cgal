// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : basic.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 

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
#else
// The following is needed for LEDA 4.4 due to min/max problems...
#  define LEDA_NO_MIN_MAX_TEMPL
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

#endif // CGAL_BASIC_H
