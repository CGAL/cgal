// Copyright (c) 2007-2008  INRIA (France).
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
//
// Author(s)     : Laurent Saboret, Pierre Alliez

#ifndef CGAL_LIGHTWEIGHT_VECTOR_3_H
#define CGAL_LIGHTWEIGHT_VECTOR_3_H

#include <CGAL/license/Poisson_surface_reconstruction_3.h>


#include <CGAL/Vector_3.h>
#include <CGAL/Origin.h>

namespace CGAL {


/// \internal
/// The Lightweight_vector_3 class represents a 3D vector (oriented).
/// The purpose of this class is to save memory as the actual vector 
/// is allocated only when needed.
///
/// \cgalModels `Kernel::Vector_3`
///
/// @tparam Gt   Geometric traits class.
template<class Gt>
class Lightweight_vector_3
{
// Public types
public:

    typedef Gt Geom_traits; ///< Geometric traits class
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::RT RT;
    typedef typename Geom_traits::Vector_3 Vector; ///< Kernel's Vector_3 class.

// Public methods
public:

    /// Vector is (0,0,0) by default.
    Lightweight_vector_3(Null_vector = NULL_VECTOR)
    {
      m_pVector = NULL;
    }
    Lightweight_vector_3(const Vector& vector)
    {
      m_pVector = new Vector(vector);
    }
    Lightweight_vector_3(FT x, FT y, FT z)
    {
      m_pVector = new Vector(x,y,z);
    }
    Lightweight_vector_3(RT hx, RT hy, RT hz, RT hw)
    {
      m_pVector = new Vector(hx,hy,hz,hw);
    }

    /// Copy constructor
    Lightweight_vector_3(const Lightweight_vector_3& that)
    {
      m_pVector = (that.m_pVector == NULL) ? NULL : new Vector(*that.m_pVector);
    }
    template <class K>
    Lightweight_vector_3(const Lightweight_vector_3<K>& that)
    {
      Vector vector = that.get_vector();
      m_pVector = (vector == NULL_VECTOR) ? NULL : new Vector(vector);
    }
    /// Operator =()
    Lightweight_vector_3& operator=(const Lightweight_vector_3& that)
    {
      if (m_pVector != NULL && that.m_pVector != NULL) 
      {
        *m_pVector = *that.m_pVector;
      }
      else
      {
        delete m_pVector;
        m_pVector = (that.m_pVector == NULL) ? NULL : new Vector(*that.m_pVector);
      }
      return *this;
    }

    /// Destructor
    ~Lightweight_vector_3()
    {
      delete m_pVector; m_pVector = NULL;
    }

    /// Compare vectors
    bool operator==(const Lightweight_vector_3& that)
    {
      return Vector(*this) == Vector(that);
    }
    bool operator!=(const Lightweight_vector_3& that)
    {
      return ! (*this == that);
    }

    /// Gets (a copy of) the actual vector. 
    operator Vector() const
    {
      if (m_pVector != NULL)
        return *m_pVector;
      else
        return NULL_VECTOR;
    }
    Vector get_vector() const
    {
      return *this;
    }

    FT x() const { return (m_pVector != NULL) ? m_pVector->x() : 0; }
    FT y() const { return (m_pVector != NULL) ? m_pVector->y() : 0; }
    FT z() const { return (m_pVector != NULL) ? m_pVector->z() : 0; }
   
    RT hx() const { return (m_pVector != NULL) ? m_pVector->hx() : 0; }
    RT hy() const { return (m_pVector != NULL) ? m_pVector->hy() : 0; }
    RT hz() const { return (m_pVector != NULL) ? m_pVector->hz() : 0; }
    RT hw() const { return (m_pVector != NULL) ? m_pVector->hw() : 1; }

    FT cartesian(int i) const
    {
      if (m_pVector != NULL)
        return m_pVector->cartesian(i);
      else
        return 0;
    }
    FT operator[](int i) const
    {
      if (m_pVector != NULL)
        return (*m_pVector)[i];
      else if (i != 3)
        return 0;
      else
        return 1;
    }
    RT homogeneous(int i) const
    {
      if (m_pVector != NULL)
        return m_pVector->homogeneous();
      else if (i != 3)
        return 0;
      else
        return 1;
    }

    int dimension() const { return 3; }

    Vector operator+(const Vector& that) const
    {
      return Vector(*this) + Vector(that);
    }
    Vector operator-(const Vector& that) const
    {
      return Vector(*this) - Vector(that);
    }
    FT operator*(const Vector& that) const
    {
      return Vector(*this) * Vector(that);
    }
    Vector operator-() const
    {
      if (m_pVector != NULL)
        return -(*m_pVector);
      else
        return NULL_VECTOR;
    }
    Vector operator/(RT c) const
    {
      if (m_pVector != NULL)
        return (*m_pVector) / c;
      else
        return NULL_VECTOR;
    }
    Vector operator*(FT c) const
    {
      if (m_pVector != NULL)
        return (*m_pVector) * c;
      else
        return NULL_VECTOR;
    }
    friend Vector operator*(FT c, const Lightweight_vector_3& vector)
    {
      return vector * c;
    }
    
    FT squared_length() const
    {
      if (m_pVector != NULL)
        return m_pVector->squared_length();
      else
        return 0;
    }

// Data
private:

    Vector* m_pVector;    // Vector (optional to save memory)
};


} //namespace CGAL

#endif //CGAL_LIGHTWEIGHT_VECTOR_3_H
