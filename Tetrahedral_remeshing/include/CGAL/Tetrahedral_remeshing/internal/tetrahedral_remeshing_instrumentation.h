// Copyright (c) 2026 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial

#ifndef CGAL_INTERNAL_TETRAHEDRAL_REMESHING_INSTRUMENTATION_H
#define CGAL_INTERNAL_TETRAHEDRAL_REMESHING_INSTRUMENTATION_H

#include <CGAL/license/Tetrahedral_remeshing.h>

// Counters used to investigate the parallel remeshing, compiled out unless
// their macro is defined. Nothing here is part of the algorithm: it exists so
// that a change to the parallel machinery can be explained, not merely timed.
//
// They are DIAGNOSTICS and never gates. Each one moves with whatever is being
// tested, and none of them predicts wall time on its own -- a zone change that
// cut retries by a quarter still ran slower, because the saving was spent
// elsewhere. Read them next to a time measurement, never instead of one.
//
//   CGAL_TR_LOCKCOUNT   lock-zone retries per elementary operation, i.e. the
//                       spin: how often a worker had to release its zone and
//                       start over because another thread held part of it.

#ifdef CGAL_TR_LOCKCOUNT
#include <cstdint>
#include <iostream>
#include <mutex>
#endif

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{

/**
* How often `apply_one()` had to give the zone back and try again -- the spin
* a lock zone pays for, counted rather than timed.
*
* Thread-local, because an atomic in this loop would add exactly the kind of
* contention it is meant to measure; the per-thread blocks are chained into a
* list at first use and summed when the program exits.
*
* A DIAGNOSTIC, not a gate: it moves with the treatment by construction, and a
* crashed run prints nothing at all, so whatever drives it must record the exit
* status separately or a dead run reads as a silent pass.
*/
struct Lockcount_counters
{
  std::uint64_t ops = 0;
  std::uint64_t retries = 0;
  Lockcount_counters* next = nullptr;
  Lockcount_counters();
};

inline std::mutex& lockcount_mutex()
{
  static std::mutex m;
  return m;
}

inline Lockcount_counters*& lockcount_head()
{
  static Lockcount_counters* head = nullptr;
  return head;
}

struct Lockcount_reporter
{
  ~Lockcount_reporter()
  {
    std::uint64_t ops = 0, retries = 0;
    std::lock_guard<std::mutex> g(lockcount_mutex());
    for (const Lockcount_counters* c = lockcount_head(); c != nullptr; c = c->next)
    {
      ops += c->ops;
      retries += c->retries;
    }
    std::cout << "LOCKCOUNT ops=" << ops << " retries=" << retries
              << " retries_per_op="
              << (ops ? double(retries) / double(ops) : 0.0) << std::endl;
  }
};

inline Lockcount_reporter& lockcount_reporter()
{
  static Lockcount_reporter r;
  return r;
}

inline Lockcount_counters::Lockcount_counters()
{
  std::lock_guard<std::mutex> g(lockcount_mutex());
  lockcount_reporter();          // must outlive every thread-local block
  next = lockcount_head();
  lockcount_head() = this;
}

inline Lockcount_counters& lockcount_counters()
{
  static thread_local Lockcount_counters c;
  return c;
}

} // end namespace internal
} // end namespace Tetrahedral_remeshing
} // end namespace CGAL

#endif // CGAL_INTERNAL_TETRAHEDRAL_REMESHING_INSTRUMENTATION_H
