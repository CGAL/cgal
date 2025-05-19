// Copyright (c) 2023 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Efi Fogel            <efifogel@gmail.com>

#ifndef CGAL_ARR_OBSERVER_H
#define CGAL_ARR_OBSERVER_H

/*! \file
 * Definition of the `Arr_observer<Arrangement>` base class mainly for backward
 * compatibility.
 */

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

namespace CGAL {

template <typename Arrangement_>
using Arr_observer = typename Arrangement_::Observer;

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
