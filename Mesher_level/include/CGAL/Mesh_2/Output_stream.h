// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESHES_OUTPUT_STREAM_H

#ifdef CGAL_MESHES_NO_OUTPUT
#  include <CGAL/IO/Verbose_ostream.h>
#  define CGAL_MESHES_OUTPUT_STREAM CGAL::Verbose_ostream()
#else
#  ifndef CGAL_MESHES_OUTPUT_ON_CERR
#    define CGAL_MESHES_OUTPUT_STREAM std::cout
#  else
#    define CGAL_MESHES_OUTPUT_STREAM std::cerr
#  endif
#endif

#endif
