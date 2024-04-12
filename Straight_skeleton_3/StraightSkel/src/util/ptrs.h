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
 * @file   util/ptrs.h
 * @author Gernot Walzl
 * @date   2012-10-30
 */

#ifndef UTIL_PTRS_H
#define UTIL_PTRS_H

#include "smarter_ptr.h"

namespace util {

class Configuration;

typedef SHARED_PTR<Configuration> ConfigurationSPtr;
typedef WEAK_PTR<Configuration> ConfigurationWPtr;

}

#endif /* UTIL_PTRS_H */
