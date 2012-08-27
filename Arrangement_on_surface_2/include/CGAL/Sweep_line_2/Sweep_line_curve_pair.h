// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_SWEEP_LINE_CURVE_PAIR_H
#define CGAL_SWEEP_LINE_CURVE_PAIR_H

/*! \file
 * Definition of the Curve_pair<Subcurve> class and related functors.
 */

namespace CGAL {

/*! \class
 * A pair of subcurves.
 */
template <class Subcurve_>
class Curve_pair
{
public:

  typedef Subcurve_                      Subcurve;

private:

  // Data members:
  std::pair<Subcurve*, Subcurve*>    m_pair;

public:

  /*! Default constructor. */
  Curve_pair(){}

  /*! Constructor from two subcurves. */
  Curve_pair (Subcurve *sc1, Subcurve *sc2)
  {
    // The smallest pointer will be the first. 
    if(sc1 < sc2)
      m_pair = std::make_pair (sc1,sc2);
    else
      m_pair = std::make_pair (sc2,sc1);
  }

  /*! Get the first subcurve. */
  Subcurve* first() const
  {
    return (m_pair.first);
  }

  /*! Get the second subcurve. */
  Subcurve* second() const
  {
    return (m_pair.second);
  }
};

/*! \struct
 * Less functor for curve pairs.
 */
template <class Subcurve_>
struct Less_curve_pair
{
  typedef Subcurve_               Subcurve;
  typedef class Curve_pair<Subcurve>   Curve_pair;

  bool operator() (const Curve_pair& pair1, const Curve_pair& pair2) const
  {
    if (pair1.first() < pair2.first())
      return true;
    if (pair1.first() > pair2.first())
      return false;
    if (pair1.second() < pair2.second())
      return true;
    return false;
  }
};

/*! \struct
 * A hash functor for curve pairs.
 */
template <class Subcurve_>
struct Curve_pair_hasher
{
  typedef Subcurve_               Subcurve;
  typedef class Curve_pair<Subcurve>   Curve_pair;

  size_t operator() (const Curve_pair& pair) const
  {
    const size_t  half_n_bits = sizeof(size_t) * 8 / 2;
    const size_t  val1 = reinterpret_cast<size_t> (pair.first());
    const size_t  val2 = reinterpret_cast<size_t> (pair.second());

    return (((val1 << half_n_bits) | (val1 >> half_n_bits)) ^ val2);
  }
};

/*! \struct
 * Equaility functor for curve pairs.
 */
template <class Subcurve_>
struct Equal_curve_pair
{
  typedef Subcurve_               Subcurve;
  typedef class Curve_pair<Subcurve>   Curve_pair;

  bool operator() (const Curve_pair& pair1, const Curve_pair& pair2) const
  {
    return (pair1.first() == pair2.first() &&
            pair1.second() == pair2.second());
  }
};

/*! \class
 * A random-access iterator that can automatically resize its container.
 */
template <class Container_>
class random_access_input_iterator
{
public:
  
  typedef Container_                                Container;
  typedef typename Container::value_type            value_type;
  typedef random_access_input_iterator<Container>   Self;

private:
  
  // Data members:
  Container          *m_container;      // The container.
  unsigned int        m_index;          // The current index.

public:

  random_access_input_iterator()
  {}

  random_access_input_iterator (Container& _container,
                                unsigned int _index = 0):
    m_container(&_container),
    m_index(_index)
  {}

  value_type& operator*()
  {
    if(m_index >= m_container->capacity())
    {
      m_container->reserve(2 * m_index + 1);
      m_container->resize(m_index+1);
    }
    else
      if(m_index >= m_container->size())
         m_container->resize(m_index+1);
    return (*m_container)[m_index];
  }

  Self& operator++()
  {
    ++m_index;
    return (*this);
  }

  Self operator++ (int)
  {
    Self temp = *this;
    ++m_index;
    return (temp);
  }

  Self& operator--()
  {
    --m_index;
    return (*this);
  }

  Self operator-- (int)
  {
    Self temp = *this;
    --m_index;
    return (temp);
  }

  bool operator== (const Self& other)
  {
    CGAL_precondition (m_container == other.m_container);
    return (m_index == other.m_index);
  }

  bool operator!=(const Self& other)
  {
    CGAL_precondition(m_container == other.m_container);
    return !(*this == other);
  }

  unsigned int operator-(const Self& other)
  {
    CGAL_precondition(m_container == other.m_container);
    return (m_index - other.m_index);
  }
};

} //namespace CGAL

#endif
