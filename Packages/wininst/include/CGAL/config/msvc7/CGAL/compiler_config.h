// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : include/CGAL/config/msvc7/CGAL/compiler_config.h
// package       : wininst
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Radu Ursu <rursu@sophia.inria.fr>
//
// ============================================================================


#if defined _MSC_VER && _MSC_VER == 1300
	#include "cl_1300.h"
#else
	#include "cl_1310.h"
#endif