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

#include "FloodFill.h"

namespace CGAL
{
namespace Qt
{

void FloodFill::operator()(
  QRgb* raw_img, uint16_t width, uint16_t x, uint16_t y, QRgb color)
{
  const QRgb old_val = raw_img[y * width + x];

  auto fill_line = [&fill_stack=(this->fill_stack), raw_img, old_val, width]
  (QRgb color, uint16_t prev_left, uint16_t prev_right, uint16_t y, int16_t dy)
  {
    QRgb* img_line = raw_img + y * width;

    for (uint16_t i = prev_left; i < prev_right; i++)
    {
      if (img_line[i] == old_val)
      {
        uint16_t j = i;

        img_line[i] = color;
        while (img_line[--j] == old_val) img_line[j] = color;
        while (img_line[++i] == old_val) img_line[i] = color;

        uint16_t left = j + 1;
        uint16_t right = static_cast<uint16_t>(i - 1);

        fill_stack.push_back({y, left, right, prev_left, prev_right, dy});
      }
    }
  };

  // first iteration
  do
  {
    uint16_t left = x;
    uint16_t right = x;

    QRgb* img_line = raw_img + y * width;
    img_line[left] = color;
    while (img_line[--left] == old_val) img_line[left] = color;
    while (img_line[++right] == old_val) img_line[right] = color;
    ++left;
    --right;

    this->fill_stack.push_back(
      {y, left, right, static_cast<uint16_t>(right + 1), right, 1});

  } while (false);

  while (!this->fill_stack.empty())
  {
    auto line = fill_stack.back();
    this->fill_stack.pop_back();

    auto ly = line.y;
    auto ll = line.left;
    auto lr = line.right;
    auto lpl = line.prev_left;
    auto lpr = line.prev_right;
    auto ldy = line.dy;

    // fill line in the new y direction
    fill_line(color, ll, lr, ly + ldy, ldy);
    // revisit the old fill line in the remaining bits
    fill_line(color, lpr + 1, lr, ly - ldy, -ldy);
    fill_line(color, ll, lpl - 1, ly - ldy, -ldy);
  }
}

} // namespace Qt
} // namespace CGAL
