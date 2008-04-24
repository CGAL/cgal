// Copyright (c) 2007  INRIA (France).
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
/// - a normal (oriented or not).
///
/// @heading Has Models:
/// - Point_with_normal_3<Geom_traits>
/// - Gyroviz_point_3<Geom_traits>

class PointWithNormal_3 : public Kernel::Point_3,
                          public DefaultConstructible, public CopyConstructible, public Assignable
{
// Public types
public:

    typedef xxx Geom_traits; ///< Kernel's geometric traits
    typedef Geom_traits::FT FT;
    typedef Geom_traits::Point_3 Point;   ///< Kernel's Point_3 class.
    typedef xxx Normal; ///< Model of OrientedNormal_3 concept.

// Public methods
public:

    /// Point is (0,0,0) by default.
    /// Normal is (0,0,0) by default.
    /// Normal is oriented by default.
    PointWithNormal_3(const Origin& o = ORIGIN);
    PointWithNormal_3(FT x, FT y, FT z);
    PointWithNormal_3(const Point& point,
                      const Normal& normal = NULL_VECTOR);

    /// Compare positions
    bool operator==(const PointWithNormal_3& that);
    bool operator!=(const PointWithNormal_3& that);

    /// Get/set normal (vector + orientation).
    const Normal& normal() const;
    Normal&       normal();
};

