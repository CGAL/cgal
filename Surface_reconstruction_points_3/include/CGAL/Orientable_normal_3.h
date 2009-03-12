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

#ifndef CGAL_ORIENTABLE_NORMAL_3_H
#define CGAL_ORIENTABLE_NORMAL_3_H

#include <CGAL/Vector_3.h>
#include <CGAL/Origin.h>

#include <boost/graph/properties.hpp>

CGAL_BEGIN_NAMESPACE


/// The Orientable_normal_3 class represents a normal vector (oriented or not).
///
/// @heading Is Model for the Concepts: Model of the OrientableNormal_3 concept.
///
/// @heading Parameters:
/// @param Gt   Kernel's geometric traits.

template<class Gt>
class Orientable_normal_3 : public Gt::Vector_3
{
// Private types
private:

  typedef typename Gt::Vector_3  Base;

// Public types
public:

    typedef Gt Geom_traits; ///< Kernel's geometric traits
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::RT RT;
    typedef typename Geom_traits::Vector_3 Vector; ///< Kernel's Vector_3 class.

// Public methods
public:

    /// Normal vector is (0,0,0) by default.
    /// Normal is oriented by default.
    Orientable_normal_3(Null_vector = NULL_VECTOR, bool oriented = true)
    : Base(NULL_VECTOR)
    {
      m_oriented = oriented;
    }
    Orientable_normal_3(const Vector& vector, bool oriented = true)
    : Base(vector)
    {
      m_oriented = oriented;
    }
    Orientable_normal_3(FT x, FT y, FT z, bool oriented = true)
    : Base(x,y,z)
    {
      m_oriented = oriented;
    }
    Orientable_normal_3(RT hx, RT hy, RT hz, RT hw, bool oriented = true)
    : Base(hx,hy,hz,hw)
    {
      m_oriented = oriented;
    }

    /// Copy constructor
    Orientable_normal_3(const Orientable_normal_3& that)
    : Base(that)
    {
      m_oriented = that.m_oriented;
    }
    template <class K>
    Orientable_normal_3(const Orientable_normal_3<K>& that)
    : Base(that)
    {
      m_oriented = that.is_oriented();
    }
    /// Operator =()
    Orientable_normal_3& operator=(const Orientable_normal_3& that)
    {
      Base::operator=(that);
      m_oriented = that.m_oriented;      
      return *this;
    }

    // Inherited operators ==() and !=() are fine.
    //bool operator==(const Orientable_normal_3& that)
    //{
    //  return ((Base&)(*this) == (Base&)that);
    //}
    //bool operator!=(const Orientable_normal_3& that)
    //{
    //  return ! (*this == that);
    //}

    /// Get (a copy of) the actual vector.
    const Vector& get_vector() const
    {
        return *this;
    }

    /// Get/set normal orientation. 
    bool is_oriented() const { return m_oriented; }
    void set_orientation(bool oriented) { m_oriented = oriented; }

// Data
private:

    bool    m_oriented; // is the normal oriented?
};


CGAL_END_NAMESPACE


namespace boost {

/// Helper type and constant to get a "vertex_normal" property map.
enum vertex_normal_t { vertex_normal } ;
BOOST_INSTALL_PROPERTY(vertex, normal);

} // namespace boost


#endif //CGAL_ORIENTABLE_NORMAL_3_H

