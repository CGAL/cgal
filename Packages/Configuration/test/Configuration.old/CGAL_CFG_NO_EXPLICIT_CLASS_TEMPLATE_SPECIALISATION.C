// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-0.9-I-05 $
// release_date  : $CGAL_Date: 1997/12/17 $
//
// file          : config/testfiles/CGAL_CFG_NO_EXPLICIT_CLASS_TEMPLATE_SPECIALISATION.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ============================================================================

// CGAL_CFG_NO_EXPLICIT_CLASS_TEMPLATE_SPECIALISATION.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler whether it knows
// the keyword typename or not.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler doesn't support specialisation of class templates,
//| the flag CGAL_CFG_NO_EXPLICIT_CLASS_TEMPLATE_SPECIALISATION is set.

template <class NT>
class Q
{
  NT d;
};

template <>
class Q<int>
{
  int d;
};

int
main()
{
  Q<double> d;
  Q<int>    i;
}

// EOF //

