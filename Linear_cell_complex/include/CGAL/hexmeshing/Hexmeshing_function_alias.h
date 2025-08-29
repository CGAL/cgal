// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>, Théo Bénard <benard320@gmail.com>
//
#ifndef HEXMESHING_FUNCTION_ALIAS_H
#define HEXMESHING_FUNCTION_ALIAS_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_outer_alias.h>
#include <functional>

namespace CGAL::internal::Hexmeshing {
  using TrimmingFunction = std::function<bool(LCC&, Dart_handle)>;
  // Identifies which 3-cell should be refined
  using MarkingFunction = std::function<bool(LCC&, Dart_handle)>;
  using DetectingFunction = std::function<bool(LCC&, Dart_handle)>;
  using DecideInsideFunction = std::function<bool(Point)>;
}



#endif