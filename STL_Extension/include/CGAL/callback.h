// Copyright (c) 2016       GeometryFactory Sarl (France)
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
// Author(s)     : Laurent Rineau
//

#ifndef CGAL_CALLBACK_H
#define CGAL_CALLBACK_H

#include <CGAL/config.h>
#  include <functional>
#  define CGAL_CALLBACK_PARAM(x) x
#  define CGAL_CALLBACK(f, ...) if(f) f(__VA_ARGS__);
#  define CGAL_CALLBACK_VAR(x) x

#endif // CGAL_CALLBACK_H
