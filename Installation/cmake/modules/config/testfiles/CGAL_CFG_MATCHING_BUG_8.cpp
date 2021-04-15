// Copyright (c) 2014  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Laurent Rineau

//| This flag is set if the compiler bugs with some "using Base::Member;" in
//| a derived class.  This bug is similar to CGAL_CFG_MATCHING_BUG_7, but
//| VC11 is not sensible to it, and VC12 says the call is ambiguous.

#define MATCHING_BUG_8 1
#include "CGAL_CFG_MATCHING_BUG_7.cpp"

// VC++12 both says that the call to `wconv(p)` is ambiguous:
//
// config\testfiles\CGAL_CFG_MATCHING_BUG_7.cpp(61) : error C2666: 'WConv::operator ()' : 2 overloads have similar conversions
//
//         config\testfiles\CGAL_CFG_MATCHING_BUG_7.cpp(51): could be 'void WConv::operator ()(const WP &)'
//
//         unable to recover from previous error(s); stopping compilation
