// Copyright (c) 2018 GeometryFactory (France).
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
// Author(s) : Simon Giraudot

#ifndef CGAL_INTERNAL_PSP_PARALLEL_CALLBACK_H
#define CGAL_INTERNAL_PSP_PARALLEL_CALLBACK_H

#include <CGAL/license/Point_set_processing_3.h>

#include <tbb/atomic.h>
#define TBB_IMPLEMENT_CPP0X 1
#include <tbb/compat/thread>

namespace CGAL {
namespace internal {
namespace Point_set_processing_3 {
  
class Parallel_callback
{
  const cpp11::function<bool(double)>& callback;
  tbb::atomic<std::size_t>& advancement;
  tbb::atomic<bool>& interrupted;
  std::size_t size;
    
public:
  Parallel_callback (const cpp11::function<bool(double)>& callback,
                     tbb::atomic<bool>& interrupted,
                     tbb::atomic<std::size_t>& advancement,
                     std::size_t size)
    : callback (callback)
    , advancement (advancement)
    , interrupted (interrupted)
    , size (size)
  { }

  void operator()()
  {
    tbb::tick_count::interval_t sleeping_time(0.00001);

    while (advancement != size)
    {
      if (!callback (advancement / double(size)))
      {
        interrupted = true;
        return;
      }
      std::this_thread::sleep_for(sleeping_time);
    }
    callback (1.);
  }
};

} // namespace Point_set_processing_3
} // namespace internal
} // namespace CGAL

#endif // CGAL_INTERNAL_PSP_PARALLEL_CALLBACK_H
