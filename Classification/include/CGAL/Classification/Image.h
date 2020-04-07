// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_IMAGE_H
#define CGAL_CLASSIFICATION_IMAGE_H

#include <CGAL/license/Classification.h>

#include <map>
#include <memory>

#define CGAL_CLASSIFICATION_IMAGE_SIZE_LIMIT 100000000

namespace CGAL {
namespace Classification {

  /// \cond SKIP_IN_MANUAL

template <typename Type>
class Image
{
  using Vector = std::vector<Type>;
  using Map = std::map<std::size_t, Type>;

  std::size_t m_width;
  std::size_t m_height;
  std::size_t m_depth;

  std::shared_ptr<Vector> m_raw;
  std::shared_ptr<Map> m_sparse;
  Type m_default;

  // Forbid using copy constructor
  Image (const Image&)
  {
  }

public:

  Image () : m_width(0), m_height(0), m_depth(0)
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
        m_raw = std::make_shared<Vector> (m_width * m_height * m_depth);
      else
        m_sparse = std::make_shared<Map> ();
    }
  }

  ~Image ()
  {
  }

  void free()
  {
    m_raw.reset();
    m_sparse.reset();
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
    if (!m_raw) // sparse case
    {
      typename Map::iterator inserted = m_sparse->insert
        (std::make_pair (coord(x,y,z), Type())).first;
      return inserted->second;
    }

    return (*m_raw)[coord(x,y,z)];
  }
  const Type& operator() (const std::size_t& x, const std::size_t& y, const std::size_t& z = 0) const
  {
    if (!m_raw) // sparse case
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
