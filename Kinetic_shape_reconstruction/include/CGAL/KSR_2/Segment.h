// Copyright (c) 2019 GeometryFactory Sarl (France).
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

#ifndef CGAL_KSR_2_SEGMENT_H
#define CGAL_KSR_2_SEGMENT_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

namespace CGAL
{

namespace KSR_2
{

template <typename GeomTraits>
class Segment
{
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;

private:

  KSR::size_t m_source;
  KSR::size_t m_target;
  KSR::size_t m_support_line;

public:

  Segment (KSR::size_t support_line) : m_support_line (support_line) { }

  const KSR::size_t& source() const { return m_source; }
  KSR::size_t& source() { return m_source; }
  const KSR::size_t& target() const { return m_target; }
  KSR::size_t& target() { return m_target; }
  KSR::size_t support_line() const { return m_support_line; }

};


}} // namespace CGAL::KSR_2


#endif // CGAL_KSR_2_POLYGON_H
