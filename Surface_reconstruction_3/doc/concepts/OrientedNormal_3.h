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


/// The OrientedNormal_3 concept represents a normal vector (oriented or not).
///
/// @heading Has Models:
/// - Oriented_normal_3<Geom_traits>

class OrientedNormal_3 : public DefaultConstructible, public CopyConstructible, public Assignable
{
// Public types
public:

    typedef xxx Vector; ///< Model of Kernel::Vector_3 concept.
    typedef typename Geom_traits::Vector_3 Vector; 

// Public methods
public:

    /// Normal vector is (0,0,0) by default.
    /// Normal is oriented by default.
    OrientedNormal_3(Null_vector = NULL_VECTOR);
    OrientedNormal_3(const Vector& vector, bool oriented = true);

    /// Get normal vector.
    Vector get_vector() const;

    /// Get normal orientation.
    bool is_oriented() const;

    /// Set normal (vector + orientation).
    void set(const Vector& vector, bool oriented = true);
};

