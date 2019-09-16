// Copyright (c) 2017 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_IMAGE_H
#define CGAL_CLASSIFICATION_IMAGE_H

#include <CGAL/license/Classification.h>

#include <boost/shared_ptr.hpp>
#include <map>

#define CGAL_CLASSIFICATION_IMAGE_SIZE_LIMIT 100000000

namespace CGAL {
namespace Classification {

  /// \cond SKIP_IN_MANUAL
  
template <typename Type>
class Image
{
  typedef std::vector<Type> Vector;
  typedef std::map<std::size_t, Type> Map;
  
  std::size_t m_width;
  std::size_t m_height;
  std::size_t m_depth;
  
  boost::shared_ptr<Vector> m_raw;
  boost::shared_ptr<Map> m_sparse;
  Type m_default;

  // Forbid using copy constructor
  Image (const Image&)
  {
  }
  
public:

  Image () : m_width(0), m_height(0), m_depth(0), m_raw (nullptr)
  {
  }
  
  Image (std::size_t width, std::size_t height, std::size_t depth = 1)
    : m_width (width)
    , m_height (height)
    , m_depth (depth)
  {
    if (m_width * m_height * m_depth > 0)
    {
      if (m_width * m_height * m_depth < CGAL_CLASSIFICATION_IMAGE_SIZE_LIMIT)
        m_raw = boost::shared_ptr<Vector> (new Vector(m_width * m_height * m_depth));
      else
        m_sparse = boost::shared_ptr<Map> (new Map());
    }
  }
  
  ~Image ()
  {
  }

  void free()
  {
    m_raw = boost::shared_ptr<Vector>();
    m_sparse = boost::shared_ptr<Map>();
  }

  Image& operator= (const Image& other)
  {
    m_raw = other.m_raw;
    m_sparse = other.m_sparse;
    m_width = other.width();
    m_height = other.height();
    m_depth = other.depth();
    return *this;
  }
  
  std::size_t width() const { return m_width; }
  std::size_t height() const { return m_height; }
  std::size_t depth() const { return m_depth; }

  inline std::size_t coord (const std::size_t& x, const std::size_t& y, const std::size_t& z) const
  {
    return z + (m_depth * y) + (m_depth * m_height * x);
  }

  Type& operator() (const std::size_t& x, const std::size_t& y, const std::size_t& z = 0)
  {
    if (m_raw == boost::shared_ptr<Vector>()) // sparse case
    {
      typename Map::iterator inserted = m_sparse->insert
        (std::make_pair (coord(x,y,z), Type())).first;
      return inserted->second;
    }

    return (*m_raw)[coord(x,y,z)];
  }
  const Type& operator() (const std::size_t& x, const std::size_t& y, const std::size_t& z = 0) const
  {
    if (m_raw == boost::shared_ptr<Vector>()) // sparse case
    {
      typename Map::iterator found = m_sparse->find (coord(x,y,z));
      if (found != m_sparse->end())
        return found->second;
      return m_default;
    }

    return (*m_raw)[coord(x,y,z)];
  }
  

};

  /// \endcond

} // namespace Classification
} // namespace CGAL



#endif // CGAL_CLASSIFICATION_IMAGE_H
