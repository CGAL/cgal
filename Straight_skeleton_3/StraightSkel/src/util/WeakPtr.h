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
 * @file   util/WeakPtr.h
 * @author Gernot Walzl
 * @date   2012-03-29
 */

#ifndef UTIL_WEAKPTR_H
#define UTIL_WEAKPTR_H

#include <memory>

namespace util {

/**
 * Adds operator== to std::weak_ptr.
 * This allows to use the WeakPtr in the STL (e.g. list).
 */
template<class T>
class WeakPtr : public std::weak_ptr<T> {
public:
    WeakPtr() : std::weak_ptr<T>() {}
    template<class Y> WeakPtr(const WeakPtr<Y>& r) : std::weak_ptr<T>(r) {}
    WeakPtr(const std::weak_ptr<T>& r) : std::weak_ptr<T>(r) {}
    WeakPtr(const std::shared_ptr<T>& r) : std::weak_ptr<T>(r) {}

    // @fixme returns true for two expired pointers, which is probably not great
    bool operator==(const std::weak_ptr<T>& other) const {
        bool result = (this->expired() && other.expired());
        if (!this->expired() && !other.expired()) {
            result = (this->lock() == other.lock());
        }
        return result;
    }
    bool operator!=(const std::weak_ptr<T>& other) const {
        return !operator==(other);
    }
};

} // namespace util

namespace std {

template<typename T>
struct owner_less<util::WeakPtr<T> >
{
    bool operator()(const util::WeakPtr<T>& lhs,
                    const util::WeakPtr<T>& rhs) const {
        return lhs.owner_before(rhs);
    }
};

} // namespace std

#endif /* UTIL_WEAKPTR_H */
