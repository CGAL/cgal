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
  boost::shared_ptr<Vector> m_raw;
  boost::shared_ptr<Map> m_sparse;
  Type m_default;

  // Forbid using copy constructor
  Image (const Image&)
  {
  }
  
public:

  Image () : m_width(0), m_height(0), m_raw (NULL)
  {
  }
  
  Image (std::size_t width, std::size_t height)
    : m_width (width),
      m_height (height)
  {
    if (m_width * m_height > 0)
    {
      if (m_width * m_height < CGAL_CLASSIFICATION_IMAGE_SIZE_LIMIT)
        m_raw = boost::shared_ptr<Vector> (new Vector(m_width * m_height));
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
    return *this;
  }
  
  std::size_t width() const { return m_width; }
  std::size_t height() const { return m_height; }

  Type& operator() (const std::size_t& x, const std::size_t& y)
  {
    if (m_raw == boost::shared_ptr<Vector>()) // sparse case
    {
      typename Map::iterator inserted = m_sparse->insert (std::make_pair (x * m_height + y, Type())).first;
      return inserted->second;
    }

    return (*m_raw)[x * m_height + y];
  }
  const Type& operator() (const std::size_t& x, const std::size_t& y) const
  {
    if (m_raw == boost::shared_ptr<Vector>()) // sparse case
    {
      typename Map::iterator found = m_sparse->find (x * m_height + y);
      if (found != m_sparse->end())
        return found->second;
      return m_default;
    }

    return (*m_raw)[x * m_height + y];
  }
  

};

  /// \endcond

} // namespace Classification
} // namespace CGAL



#endif // CGAL_CLASSIFICATION_IMAGE_H
