// Copyright (c) 2017 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_ATTRIBUTE_EFFECT_H
#define CGAL_CLASSIFICATION_ATTRIBUTE_EFFECT_H

namespace CGAL {

namespace Classification {

namespace Attribute {

enum Effect /// Defines the effect of an attribute on a type.
    {
      FAVORING = 0, ///< High values of the attribute favor this type
      NEUTRAL = 1, ///< The attribute has no effect on this type
      PENALIZING = 2 ///< Low values of the attribute favor this type
    };
  
} // Attribute

} // Classification

} // CGAL

#endif
