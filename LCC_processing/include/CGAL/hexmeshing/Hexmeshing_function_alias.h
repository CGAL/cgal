// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>, Théo Bénard <benard320@gmail.com>
//

#ifndef CGAL_HEXMESHING_FUNCTION_ALIAS_H
#define CGAL_HEXMESHING_FUNCTION_ALIAS_H

#include <CGAL/license/LCC_processing.h>

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <functional>

namespace CGAL::internal::Hexmeshing
{
  using TrimmingFunction = std::function<bool(LCC&, Dart_descriptor)>;
  // Identifies which 3-cell should be refined
  using MarkingFunction = std::function<bool(LCC&, Dart_descriptor)>;
  using DetectingFunction = std::function<bool(LCC&, Dart_descriptor)>;
  using DecideInsideFunction = std::function<bool(Point)>;
}

#endif
