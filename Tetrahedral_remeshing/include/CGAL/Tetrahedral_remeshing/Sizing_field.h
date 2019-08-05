// Copyright (c) 2019 GeometryFactory (France).
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
// Author(s)     : Jane Tournois

#ifndef CGAL_SIZING_FIELD_H
#define CGAL_SIZING_FIELD_H

namespace CGAL
{
  /*!
  * Sizing field virtual class
  */
  template <class Kernel>
  class Sizing_field
  {
  public:
    typedef Kernel                      K;
    typedef typename Kernel::FT         FT;
    typedef typename Kernel::Point_3    Point_3;

  public:
    virtual FT operator()(const Point_3& p) const = 0;
  };

}//end namespace CGAL

#endif //CGAL_SIZING_FIELD_H
