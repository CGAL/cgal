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

#include <CGAL/function.h>

#include <CGAL/thread.h>

namespace CGAL {
namespace internal {
namespace Point_set_processing_3 {
  
class Parallel_callback
{
  const cpp11::function<bool(double)>& m_callback;
  cpp11::atomic<std::size_t>* m_advancement;
  cpp11::atomic<bool>* m_interrupted;
  std::size_t m_size;
  bool m_creator;
  cpp11::thread* m_thread;

  // assignment operator shouldn't be used (m_callback is const ref)
  Parallel_callback& operator= (const Parallel_callback&)
  {
    return *this;
  }
  
public:
  Parallel_callback (const cpp11::function<bool(double)>& callback,
                     std::size_t size,
                     std::size_t advancement = 0,
                     bool interrupted = false)
    : m_callback (callback)
    , m_advancement (new cpp11::atomic<std::size_t>())
    , m_interrupted (new cpp11::atomic<bool>())
    , m_size (size)
    , m_creator (true)
    , m_thread (NULL)
  {
    // cpp11::atomic only has default constructor, initialization done in two steps
    *m_advancement = advancement;
    *m_interrupted = interrupted;
    if (m_callback)
      m_thread = new cpp11::thread (*this);
  }

  Parallel_callback (const Parallel_callback& other)
    : m_callback (other.m_callback)
    , m_advancement (other.m_advancement)
    , m_interrupted (other.m_interrupted)
    , m_size (other.m_size)
    , m_creator (false)
    , m_thread (NULL)
  {

  }


  ~Parallel_callback ()
  {
    if (m_creator)
    {
      delete m_advancement;
      delete m_interrupted;
    }
    if (m_thread != NULL)
      delete m_thread;
  }

  cpp11::atomic<std::size_t>& advancement() { return *m_advancement; }
  cpp11::atomic<bool>& interrupted() { return *m_interrupted; }
  void join()
  {
    if (m_thread != NULL)
      m_thread->join();
  }

  void operator()()
  {
    while (*m_advancement != m_size)
    {
      if (!m_callback (*m_advancement / double(m_size)))
        *m_interrupted = true;
      if (*m_interrupted)
        return;
      cpp11::sleep_for (0.00001);
    }
    m_callback (1.);
  }
};

} // namespace Point_set_processing_3
} // namespace internal
} // namespace CGAL

#endif // CGAL_INTERNAL_PSP_PARALLEL_CALLBACK_H
