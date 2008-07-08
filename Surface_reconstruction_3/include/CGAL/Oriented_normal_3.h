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

#ifndef CGAL_ORIENTED_NORMAL_3_H
#define CGAL_ORIENTED_NORMAL_3_H

#include <CGAL/Vector_3.h>
#include <CGAL/Origin.h>

#include <boost/graph/properties.hpp>

CGAL_BEGIN_NAMESPACE


/// The Oriented_normal_3 class represents a normal vector (oriented or not).
/// The normal vector is allocated only when needed.
///
/// @heading Is Model for the Concepts: Model of the OrientedNormal_3 concept.
///
/// @heading Parameters:
/// @param Gt   Kernel's geometric traits.

template<class Gt>
class Oriented_normal_3
{
// Public types
public:

    typedef Gt Geom_traits; ///< Kernel's geometric traits
    typedef typename Geom_traits::Vector_3 Vector; 

// Public methods
public:

    /// Normal vector is (0,0,0) by default.
    /// Normal is oriented by default.
    Oriented_normal_3(Null_vector = NULL_VECTOR)
    {
      m_pNormal = NULL;
      m_oriented = true;
    }
    Oriented_normal_3(const Vector& vector, bool oriented = true)
    {
      m_pNormal = new Vector(vector);
      m_oriented = oriented;
    }

    /// Copy constructor
    Oriented_normal_3(const Oriented_normal_3& that)
    {
      m_pNormal = (that.m_pNormal == NULL) ? NULL : new Vector(*that.m_pNormal);
      m_oriented = that.m_oriented;
    }
    /// Operator =()
    Oriented_normal_3& operator=(const Oriented_normal_3& that)
    {
      if (m_pNormal != NULL && that.m_pNormal != NULL) 
      {
        *m_pNormal = *that.m_pNormal;
      }
      else
      {
        delete m_pNormal;
        m_pNormal = (that.m_pNormal == NULL) ? NULL : new Vector(*that.m_pNormal);
      }

      m_oriented = that.m_oriented;
      
      return *this;
    }

    /// Destructor
    ~Oriented_normal_3()
    {
      delete m_pNormal;
    }

    /// Get normal vector. 
    Vector get_vector() const
    {
      if(m_pNormal != NULL)
        return *m_pNormal;
      else
        return CGAL::NULL_VECTOR;
    }

    /// Get normal orientation. 
    bool is_oriented() const { return m_oriented; }

    /// Set normal (vector + orientation). 
    void set(const Vector& vector, bool oriented = true)
    {
      if(m_pNormal == NULL)
        m_pNormal = new Vector(vector);
      else
        *m_pNormal = vector;
      m_oriented = oriented;
    }

// Data
private:

    // PA: why is normal optional here? to save memory?
    Vector *m_pNormal;    // normal vector (optional)
    bool    m_oriented;
};


CGAL_END_NAMESPACE


namespace boost {

/// Helper type and constant to get a "vertex_normal" property map.
enum vertex_normal_t { vertex_normal } ;
BOOST_INSTALL_PROPERTY(vertex, normal);

} // namespace boost


#endif //CGAL_ORIENTED_NORMAL_3_H

