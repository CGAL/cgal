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
#ifndef HEXMESHING_GENERIC_POINT_H
#define HEXMESHING_GENERIC_POINT_H

#include <array>


namespace CGAL::internal::Hexmeshing {
  /**
   * @brief A generic 3D point class with basic arithmetic operations
   * @tparam T The numeric type used for coordinates (e.g., int, char)
   *
   * This class represents a point in 3D space with x, y, and z coordinates.
   * It provides basic arithmetic operations such as addition, subtraction,
   * and division. The class is designed to be used with different numeric
   * types through template parameter T.
   */
  template <typename T>
  struct GenericPointForHexmeshing {
    T x, y, z;  ///< The x, y, z coordinates of the point

    /// @brief Default constructor, initializes coordinates to 0
    constexpr GenericPointForHexmeshing(): x(0), y(0), z(0) {}

    /// @brief Constructor with explicit coordinates
    /// @param x The x coordinate
    /// @param y The y coordinate
    /// @param z The z coordinate
    constexpr GenericPointForHexmeshing(T x, T y, T z): x(x), y(y), z(z) {}

    /// @brief Constructor from array of 3 integers
    /// @param l Array containing the coordinates [x,y,z]
    constexpr GenericPointForHexmeshing(std::array<int,3> l): x(l[0]), y(l[1]), z(l[2]) {}
    GenericPointForHexmeshing operator+(const GenericPointForHexmeshing& p) const { return { x + p.x, y + p.y, z + p.z }; }
    GenericPointForHexmeshing operator/(const GenericPointForHexmeshing& p) const { return { x / p.x, y / p.y, z / p.z }; }
    GenericPointForHexmeshing operator-(const GenericPointForHexmeshing& p) const { return { x - p.x, y - p.y, z - p.z }; }
    GenericPointForHexmeshing operator-() const { return {static_cast<T>(-x), static_cast<T>(-y), static_cast<T>(-z)}; }
    GenericPointForHexmeshing& operator+=(const GenericPointForHexmeshing& p) { *this = *this + p; return *this; }
    GenericPointForHexmeshing& operator-=(const GenericPointForHexmeshing& p) { *this = *this - p; return *this; }
    GenericPointForHexmeshing& operator/=(const GenericPointForHexmeshing& p) { *this = *this / p; return *this; }
    T& operator[](int i) { assert(i >= 0 && i <= 2); return i == 0 ? x : i == 1 ? y : i == 2 ? z : z; };
    const T& operator[](int i) const { assert(i >= 0 && i <= 2); return i == 0 ? x : i == 1 ? y : i == 2 ? z : z; };
    bool operator==(const GenericPointForHexmeshing& p) const { return p.x == x && p.y == y && p.z == z; }
    bool operator!=(const GenericPointForHexmeshing& p) const { return p.x != x or p.y != y or p.z != z; }

    /// @brief Type conversion operator
    /// @tparam K Target numeric type
    /// @return Point with coordinates converted to type K
    template <typename K>
    operator GenericPointForHexmeshing<K>() const { return {static_cast<K>(x), static_cast<K>(y), static_cast<K>(z)}; }
  };

  using PointInt = GenericPointForHexmeshing<int>;
  using PointChar = GenericPointForHexmeshing<char>;
  using AreaId = PointChar;
}


#endif