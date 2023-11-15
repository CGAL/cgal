// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSP_2_SEGMENT_H
#define CGAL_KSP_2_SEGMENT_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

namespace CGAL
{

namespace KSP_2
{

class Segment
{
private:

  std::size_t m_input_idx;
  std::size_t m_source_idx;
  std::size_t m_target_idx;
  std::size_t m_support_line_idx;

public:

  Segment () { }

  Segment (std::size_t input_idx, std::size_t support_line_idx)
    : m_input_idx (input_idx), m_support_line_idx (support_line_idx) { }

  const std::size_t& input_idx() const { return m_input_idx; }
  std::size_t& input_idx() { return m_input_idx; }
  const std::size_t& source_idx() const { return m_source_idx; }
  std::size_t& source_idx() { return m_source_idx; }
  const std::size_t& target_idx() const { return m_target_idx; }
  std::size_t& target_idx() { return m_target_idx; }
  const std::size_t& support_line_idx() const { return m_support_line_idx; }
  std::size_t& support_line_idx() { return m_support_line_idx; }
};


}} // namespace CGAL::KSP_2


#endif // CGAL_KSP_2_POLYGON_H
