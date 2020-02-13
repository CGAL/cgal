// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_PSP_INTERNAL_PARALLEL_CALLBACK_H
#define CGAL_PSP_INTERNAL_PARALLEL_CALLBACK_H

#include <CGAL/license/Point_set_processing_3.h>

#include <functional>

#include <CGAL/thread.h>

namespace CGAL {
namespace Point_set_processing_3 {
namespace internal {
  
class Parallel_callback
{
  const std::function<bool(double)>& m_callback;
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
  Parallel_callback (const std::function<bool(double)>& callback,
                     std::size_t size,
                     std::size_t advancement = 0,
                     bool interrupted = false)
    : m_callback (callback)
    , m_advancement (new cpp11::atomic<std::size_t>())
    , m_interrupted (new cpp11::atomic<bool>())
    , m_size (size)
    , m_creator (true)
    , m_thread (nullptr)
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
    , m_thread (nullptr)
  {

  }


  ~Parallel_callback ()
  {
    if (m_creator)
    {
      delete m_advancement;
      delete m_interrupted;
    }
    if (m_thread != nullptr)
      delete m_thread;
  }

  cpp11::atomic<std::size_t>& advancement() { return *m_advancement; }
  cpp11::atomic<bool>& interrupted() { return *m_interrupted; }
  void join()
  {
    if (m_thread != nullptr)
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

} // namespace internal
} // namespace Point_set_processing_3
} // namespace CGAL

#endif // CGAL_PSP_INTERNAL_PARALLEL_CALLBACK_H
