// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
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
// file          : config/testfiles/CGAL_CFG_NOMINMAX_BUG.C
// package       : Configuration (2.3)
// author(s)     : Radu Ursu
//
// coordinator   : Radu Ursu -- rursu@sophia.inria.fr
//
// ======================================================================

// CGAL_CFG_NOMINMAX_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This is a test-case for a bug in VC++ 7.0 that redefines min(a, b) and max(a, b)
//| Files concerned: windows.h, windef.h
//| When the bug is present, CGAL_CFG_NOMINMAX_BUG is set
//| The file basic.h should check if this bug is present and if so, define NOMINMAX flag

#if defined _MSC_VER
#error "NOMINMAX flag should be set"
#endif
