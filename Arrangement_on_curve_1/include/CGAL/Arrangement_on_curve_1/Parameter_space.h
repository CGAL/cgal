// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_PARAMETER_SPACE_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_PARAMETER_SPACE_H

namespace CGAL {
namespace Arrangement_on_curve_1 {

enum class Interval_parameter_space {
  LEFT_BOUNDARY = 0,
  RIGHT_BOUNDARY,
  INTERIOR
};

using Parameter_space = Interval_parameter_space;

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
