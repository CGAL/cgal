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
 * @file   smarter_ptr.h
 * @author Gernot Walzl
 * @date   2012-10-17
 */

#ifndef SMARTER_PTR_H
#define SMARTER_PTR_H

#include "config.h"

#include <algorithm>
#include <iterator>
#include <memory>

# define SHARED_PTR std::shared_ptr
# define WEAK_PTR std::weak_ptr

namespace util {

template<class Iter, class V>
auto weak_find(Iter begin, Iter end, const std::weak_ptr<V>& target)
{
    return std::find_if(begin, end, [&](const std::weak_ptr<V>& w) {
                                        auto sp = w.lock();
                                        auto tp = target.lock();
                                        return sp && tp && (sp == tp);
                                      });
}

} // namespace util

#endif /* SMARTER_PTR_H */
