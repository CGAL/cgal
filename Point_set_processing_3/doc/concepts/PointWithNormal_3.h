// Copyright (c) 2007-2008  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Laurent Saboret, Pierre Alliez


/// The PointWithNormal_3 concept represents a 3D point with:
/// - a position
/// - a normal (oriented).
///
/// @heading Has Models: Point_with_normal_3<GeomTraits>

class PointWithNormal_3 
  : public Kernel::Point_3,
    public DefaultConstructible, public CopyConstructible, public Assignable, public EqualityComparable
{
// Public types
public:

    typedef xxx Geom_traits; ///< Kernel's geometric traits
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::RT RT;
    typedef typename Geom_traits::Point_3  Point;  ///< Model of Kernel::Point_3.
    typedef typename Geom_traits::Vector_3 Vector; ///< Model of Kernel::Vector_3.

// Public methods
public:

    /// Point is (0,0,0) by default.
    /// Normal is (0,0,0) by default.
    PointWithNormal_3(const Origin& o = ORIGIN);
    PointWithNormal_3(FT x, FT y, FT z,
                      const Vector& normal = NULL_VECTOR);
    PointWithNormal_3(RT hx, RT hy, RT hz, RT hw,
                      const Vector& normal = NULL_VECTOR);
    PointWithNormal_3(const Point& point,
                      const Vector& normal = NULL_VECTOR);

    /// Get/set position.
    const Point& position() const;
    Point&       position();

    /// Get/set normal.
    const Vector& normal() const;
    Vector&       normal();
};

