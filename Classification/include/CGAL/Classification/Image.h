// Copyright (c) 2016  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_IMAGE_H
#define CGAL_CLASSIFICATION_IMAGE_H

namespace CGAL {
namespace Classification {

  /// \cond SKIP_IN_MANUAL
  
template <typename Type>
class Image
{
  std::size_t m_width;
  std::size_t m_height;
  Type* m_raw;


public:

  Image () : m_width(0), m_height(0), m_raw (NULL)
  {
  }
  
  Image (std::size_t width, std::size_t height)
    : m_width (width),
      m_height (height)
  {
    if (m_width * m_height > 0)
      m_raw = new Type[width * height]();
    else
      m_raw = NULL;
  }
  
  ~Image ()
  {
    free();
  }

  void free()
  {
    if (m_raw != NULL)
      delete[] m_raw;
    m_raw = NULL;
  }

  Image (const Image& other)
    : m_width (other.width()),
      m_height (other.height())

  {
    if (m_width * m_height > 0)
      {
        m_raw = new Type[m_width * m_height];
        std::copy (other.m_raw, other.m_raw + (m_width * m_height), this->m_raw);
      }
    else
      m_raw = NULL;
  }
  Image& operator= (const Image& other)
  {
    if (m_raw != NULL)
      delete[] m_raw;

    m_raw = NULL;
    m_width = other.width();
    m_height = other.height();
    if (m_width * m_height > 0)
      {
        m_raw = new Type[m_width * m_height];
        std::copy (other.m_raw, other.m_raw + (m_width * m_height), this->m_raw);
      }
    
    return *this;
  }
  
  std::size_t width() const { return m_width; }
  std::size_t height() const { return m_height; }

  Type& operator() (const std::size_t& x, const std::size_t& y)
  {
    //    return m_raw[y * m_width + x];
    return m_raw[x * m_height + y];
  }
  const Type& operator() (const std::size_t& x, const std::size_t& y) const
  {
    //    return m_raw[y * m_width + x];
    return m_raw[x * m_height + y];
  }
  

};

  /// \endcond

} // namespace Classification
} // namespace CGAL



#endif // CGAL_CLASSIFICATION_IMAGE_H
