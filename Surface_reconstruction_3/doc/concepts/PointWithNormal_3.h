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


/// A PointWithNormal_3 concept represents a 3D point with:
/// - a position
/// - a normal (oriented).
///
/// @heading Has Models:
/// - Point_with_normal_3<Geom_traits, Normal_3>

class PointWithNormal_3 
  : public Kernel::Point_3,
    public DefaultConstructible, public CopyConstructible, public Assignable, public EqualityComparable
{
// Public types
public:

    typedef xxx Geom_traits; ///< Kernel's geometric traits
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::RT RT;
    typedef typename Geom_traits::Point_3  Point;  ///< Kernel's Point_3 class.
    typedef Normal_3 Normal; ///< Model of Kernel::Vector_3 or of OrientableNormal_3.

// Public methods
public:

    /// Point is (0,0,0) by default.
    /// Normal is (0,0,0) by default.
    PointWithNormal_3(const Origin& o = ORIGIN);
    PointWithNormal_3(FT x, FT y, FT z,
                      const Normal& normal = NULL_VECTOR);
    PointWithNormal_3(RT hx, RT hy, RT hz, RT hw,
                      const Normal& normal = NULL_VECTOR);
    PointWithNormal_3(const Point& point,
                      const Normal& normal = NULL_VECTOR);

    /// Set position.
    void set_position(const Point& point);

    /// Get/set normal.
    const Normal& normal() const;
    Normal&       normal();
};

