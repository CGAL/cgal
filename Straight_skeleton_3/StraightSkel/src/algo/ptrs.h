// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   algo/ptrs.h
 * @author Gernot Walzl
 * @date   2012-02-06
 */

#ifndef ALGO_PTRS_H
#define ALGO_PTRS_H

#include "smarter_ptr.h"

namespace algo {

class Controller;

typedef SHARED_PTR<Controller> ControllerSPtr;
typedef WEAK_PTR<Controller> ControllerWPtr;

}

#endif /* ALGO_PTRS_H */
