// Copyright (c) 2007,2009,2011 Max-Planck-Institute for Computer Science (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>

#ifndef CGAL_ARR_OVERLAY_H
#define CGAL_ARR_OVERLAY_H

/*! \file
 * Helping file to include Arr_overlay_2 for backward compatibility.
 */

#if (defined __GNUC__)
  #warning Arr_overlay.h is DEPRECATED, please include Arr_overlay_2.h instead
#elif (defined _MSC_VER)
  #pragma message("Arr_overlay.h is DEPRECATED, please include Arr_overlay_2.h instead")
#endif

#include <CGAL/Arr_overlay_2.h>

#endif
