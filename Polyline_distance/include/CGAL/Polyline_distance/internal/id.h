// Copyright (c) 2019 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//
// =============================================================================

#pragma once

#include <cstdint>
#include <functional>
#include <limits>


namespace CGAL {
namespace Polyline_distance {
namespace internal {
// Typesafe ID class such that there are compiler errors if different IDs are
// mixed. The template parameter T is just there to assure this behavior.
// Additionally, we have a member function which can check for validity.
template <typename T>
struct ID {
public:
    using IDType = uint32_t;
    static constexpr IDType invalid_value = std::numeric_limits<IDType>::max();

    ID(IDType id = invalid_value) : id(id) {}

    operator IDType() const { return id; }
    IDType operator+(ID<T> other) const { return id + other.id; }
    IDType operator+(int offset) const { return id + offset; }
    IDType operator+(size_t offset) const { return id + offset; }
    IDType operator-(ID<T> other) const { return id - other.id; }
    IDType operator-(int offset) const { return id - offset; }
    IDType operator/(int div) const { return id / div; }
    IDType operator+=(ID<T> other) { return id += other.id; }
    IDType operator-=(ID<T> other) { return id -= other.id; }
    IDType operator=(ID<T> other) { return id = other.id; }
    IDType operator++() { return ++id; }
    IDType operator--() { return --id; }
    // FIXME:
    // bool operator==(ID<T> other) const { return id == other.id; }
    // bool operator==(IDType other) const { return id == other; }
    bool operator!=(ID<T> other) const { return id != other.id; }

    bool valid() const { return id != invalid_value; }
    void invalidate() { id = invalid_value; }

private:
    IDType id;
};

} // namespace internal
} // namespace Polyline_distance
} // namespace CGAL

// define custom hash function to be able to use IDs with maps/sets
namespace std
{

template <typename T>
struct hash<CGAL::Polyline_distance::internal::ID<T>> {
    using IDType = typename CGAL::Polyline_distance::internal::ID<T>::IDType;
    std::size_t operator()(CGAL::Polyline_distance::internal::ID<T> const& id) const noexcept
    {
        return std::hash<IDType>()(id);
    }
};

}  // namespace std
