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


/// The OrientableNormal_3 concept represents a normal vector (oriented or not).
///
/// @heading Has Models:
/// - Orientable_normal_3<Geom_traits>

class OrientableNormal_3
  : public Kernel::Vector_3,
    public DefaultConstructible, public CopyConstructible, public Assignable, public EqualityComparable
{
// Public types
public:

    typedef xxx Geom_traits; ///< Kernel's geometric traits
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::RT RT;
    typedef typename Geom_traits::Vector_3 Vector; ///< Kernel's Vector_3 class.

// Public methods
public:

    /// Normal vector is (0,0,0) by default.
    /// Normal is oriented by default.
    OrientableNormal_3(Null_vector = NULL_VECTOR);
    OrientableNormal_3(const Vector& vector, bool oriented = true);
    OrientableNormal_3(FT x, FT y, FT z, bool oriented = true);
    OrientableNormal_3(RT hx, RT hy, RT hz, RT hw, bool oriented = true);

    /// Get (a copy of) the actual vector.
    Vector get_vector() const;

    /// Get/set normal orientation. 
    bool is_oriented() const;
    void set_orientation(bool oriented);
};

