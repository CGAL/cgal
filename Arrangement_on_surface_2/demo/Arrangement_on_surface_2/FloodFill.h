// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ahmed Essam <theartful.ae@gmail.com>

#ifndef ARRANGEMENT_DEMO_FLOOD_FILL_H
#define ARRANGEMENT_DEMO_FLOOD_FILL_H

#include <QColor>
#include <vector>

namespace CGAL
{
namespace Qt
{

// scanline flood fill
class FloodFill
{
public:
  // this currently assumes that there is a "border" in the boundaries that
  // will prevent the flood fill from going there
  // this way we don't check bounadry conditions!
  void
  operator()(QRgb* raw_img, uint16_t width, uint16_t x, uint16_t y, QRgb color);

private:
  struct FillLine
  {
    uint16_t y;
    uint16_t left;
    uint16_t right;
    uint16_t prev_left;
    uint16_t prev_right;
    // y = prev_y + dy
    int16_t dy;

    FillLine(
      uint16_t y_, uint16_t left_, uint16_t right_, uint16_t prev_left_,
      uint16_t prev_right_, int16_t dy_) :
        y{y_},
        left{left_}, right{right_}, prev_left{prev_left_},
        prev_right{prev_right_}, dy{dy_}
    {
    }
  };

  std::vector<FillLine> fill_stack;
};

} // namespace Qt
} // namespace CGAL

#endif
