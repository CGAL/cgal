// Copyright (c) 2026
// GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Andreas Fabri

#ifndef CGAL_MONOSTATE_H
#define CGAL_MONOSTATE_H

namespace CGAL {

class monostate {
  void* ptr = nullptr;
};

constexpr bool operator==(monostate, monostate) noexcept
{
  return true;
}

constexpr bool operator!=(monostate, monostate) noexcept
{
  return false;
}

constexpr bool operator<(monostate, monostate) noexcept
{
  return false;
}

constexpr bool operator>(monostate, monostate) noexcept
{
  return false;
}

constexpr bool operator<=(monostate, monostate) noexcept
{
  return true;
}

constexpr bool operator>=(monostate, monostate) noexcept
{
  return true;
}

} // namespace CGAL

#endif