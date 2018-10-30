// Copyright (c) 2018  Carnegie Mellon University (USA), GeometryFactory (France)
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s) : Christina Vaz, Keenan Crane, Andreas Fabri

#ifndef CGAL_HEAT_METHOD_3_INTERNAL_V2V_H
#define CGAL_HEAT_METHOD_3_INTERNAL_V2V_H

#include <CGAL/license/Heat_method_3.h>

#ifndef DOXYGEN_RUNNING

namespace CGAL {
namespace Heat_method_3 {

  template <typename TM>
  struct V2V {

    V2V(const TM&)
    {}

    template <typename T>
    const T& operator()(const T& t) const
    {
      return t;
    }
  };
#endif
} // namespace Heat_method_3
} // namespace CGAL
#endif //CGAL_HEAT_METHOD_3_INTERNAL_V2V_H
