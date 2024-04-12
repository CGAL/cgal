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
 * @file   util/Timer.h
 * @author Gernot Walzl
 * @date   2013-06-25
 */

#ifndef UTIL_TIMER_H
#define UTIL_TIMER_H

#include <boost/date_time/posix_time/posix_time.hpp>

namespace util {

class Timer {
public:
    virtual ~Timer();
    static double now();
protected:
    Timer();
};

}

#endif /* UTIL_TIMER_H */

