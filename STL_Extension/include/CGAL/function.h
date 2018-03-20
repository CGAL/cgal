// Copyright (c) 2017 GeometryFactory Sarl (France).
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_FUNCTION_H
#define CGAL_FUNCTION_H

#include <CGAL/config.h>
#ifndef CGAL_CFG_NO_STD_FUNCTION
#  include <functional>
#else
#  include <boost/function.hpp>
#endif

namespace CGAL {

namespace cpp11 {

#ifndef CGAL_CFG_NO_STD_FUNCTION
using std::function;
#else
using boost::function;
#endif

} // cpp11

} //namespace CGAL

#endif // CGAL_FUNCTION_H
