// Copyright (c) 2014  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
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
