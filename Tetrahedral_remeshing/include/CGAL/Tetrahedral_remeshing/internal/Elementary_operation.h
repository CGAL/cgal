// Copyright (c) 2025 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Iasonas Manolas, Jane Tournois

#ifndef CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H
#define CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_instrumentation.h>

#include <CGAL/tags.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/blocked_range.h>
#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/task_arena.h>
#endif

#include <algorithm>
#include <atomic>
#include <iterator>
#include <random>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#include <CGAL/Real_timer.h>
#include <cstddef>
#include <iostream>
#endif

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {


template <typename C3t3_, typename ElementType, typename ElementRange>
class Elementary_operation
{
public:
  using C3t3 = C3t3_;
  using Triangulation = typename C3t3::Triangulation;
  using Element_type = ElementType;
  using Element_range = ElementRange;

  Elementary_operation() = default;
  virtual ~Elementary_operation() = default;

  virtual Element_range get_elements(const C3t3& c3t3) const = 0;
  virtual bool execute_operation(const Element_type& e, C3t3& c3t3) = 0;
  virtual std::string operation_name() const = 0;
};

template <typename Operation>
class Elementary_operation_execution_sequential
{
public:
  using C3t3 = typename Operation::C3t3;
  using Element_range = typename Operation::Element_range;

  bool execute(Operation& op, C3t3& c3t3) const
  {
    Element_range candidates = op.get_elements(c3t3);
    if (candidates.empty())
      return false;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::size_t nb_done = 0;
    const std::size_t nb_candidates = candidates.size();
    CGAL::Real_timer timer;
    timer.start();
#endif
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
    std::size_t nb_processed = 0;
#endif
    for (const auto& element : candidates)
    {
      if (op.execute_operation(element, c3t3))
      {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        ++nb_done;
#endif
      }
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
      std::cout << "\r" << op.operation_name() << "... ("
                << ++nb_processed << "/" << nb_candidates << ")";
      std::cout.flush();
#endif
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    timer.stop();
    std::cout << op.operation_name() << ": " << nb_done << "/"
              << nb_candidates << " done ("
              << timer.time() << " sec)." << std::endl;
#endif
    return true;
  }
};

/**
* Picks a container type from the triangulation's concurrency tag:
* `SequentialContainer` when remeshing sequentially, `ConcurrentContainer`
* when remeshing in parallel. The sequential type is spelled out at every use
* site, so that making an operation parallelizable cannot silently change the
* container the sequential path uses.
*/
template <typename ConcurrencyTag,
          typename SequentialContainer,
          typename ConcurrentContainer>
struct Concurrency_selected_container
{
  using type = SequentialContainer;
};

#ifdef CGAL_LINKED_WITH_TBB
template <typename SequentialContainer, typename ConcurrentContainer>
struct Concurrency_selected_container<CGAL::Parallel_tag,
                                      SequentialContainer,
                                      ConcurrentContainer>
{
  using type = ConcurrentContainer;
};
#endif

template <typename ConcurrencyTag,
          typename SequentialContainer,
          typename ConcurrentContainer>
using Concurrency_selected_container_t =
  typename Concurrency_selected_container<ConcurrencyTag,
                                          SequentialContainer,
                                          ConcurrentContainer>::type;

/**
* The parallel counterpart of `Elementary_operation_execution_sequential`.
*
* An operation is parallelizable when, in addition to the `Elementary_operation`
* interface, it provides
*   - `bool lock_zone(const Element_type&, const C3t3&) const`, which locks
*     every element `execute_operation()` may touch, and
*   - `static constexpr bool requires_ordered_processing`, which says whether
*     the order of `get_elements()` carries meaning.
*
* Neither is a virtual of `Elementary_operation`: this class is a template, so
* they are found on the concrete operation. The sequential path never names
* them, and is therefore unaffected by parallelism being available.
*
* `requires_ordered_processing == true` drains a concurrent queue in the order
* `get_elements()` produced, so that threads still take the most-wanted
* elements first. `false` shuffles instead, to spread the threads over the
* triangulation and keep lock conflicts down.
*/
#ifdef CGAL_LINKED_WITH_TBB
template <typename Operation>
class Elementary_operation_execution_parallel
{
public:
  using C3t3 = typename Operation::C3t3;
  using Element_type = typename Operation::Element_type;
  using Element_range = typename Operation::Element_range;

  bool execute(Operation& op, C3t3& c3t3) const
  {
    std::vector<Element_type> candidates = collect(op, c3t3);
    if (candidates.empty())
      return false;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    CGAL::Real_timer timer;
    timer.start();
    const std::size_t nb_candidates = candidates.size();
#endif

    if constexpr (Operation::requires_ordered_processing)
      run_ordered(candidates, op, c3t3);
    else
      run_unordered(candidates, op, c3t3);

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    timer.stop();
    std::cout << op.operation_name() << ": " << nb_candidates
              << " candidates (" << timer.time() << " sec, parallel)."
              << std::endl;
#endif
    return true;
  }

private:
  static std::vector<Element_type> collect(const Operation& op, const C3t3& c3t3)
  {
    Element_range range = op.get_elements(c3t3);
    if constexpr (std::is_same_v<Element_range, std::vector<Element_type>>)
      return std::move(range);
    else
    {
      std::vector<Element_type> candidates;
      candidates.reserve(std::distance(range.begin(), range.end()));
      std::copy(range.begin(), range.end(), std::back_inserter(candidates));
      return candidates;
    }
  }

  /**
  * Locks the zone of `element`, retrying until it succeeds, then runs the
  * operation and releases the zone. The retry is unbounded on purpose: a
  * failed lock means another thread holds part of this zone, and every zone
  * is released as soon as its operation ends, so the wait is finite.
  */
  static void apply_one(const Element_type& element, Operation& op, C3t3& c3t3)
  {
#ifdef CGAL_TR_LOCKCOUNT
    Lockcount_counters& lc = lockcount_counters();
    ++lc.ops;
#endif
    while (!op.lock_zone(element, c3t3))
    {
#ifdef CGAL_TR_LOCKCOUNT
      ++lc.retries;
#endif
      c3t3.triangulation().unlock_all_elements();
      std::this_thread::yield();
    }
    op.execute_operation(element, c3t3);
    c3t3.triangulation().unlock_all_elements();
  }

  /**
  * How many elements a worker of an ORDERED operation may hold back before it
  * stops taking new ones.
  *
  * An ordered operation cannot defer to the end of the pass: its candidate
  * order carries meaning -- shortest edge first for the collapse, longest for
  * the split -- and postponing an element until everything else is done
  * changes the trajectory rather than the schedule. What it can do is stop
  * WAITING: a worker that cannot take a zone sets the element aside, takes the
  * next candidate, and comes back to the set-aside ones as soon as a few have
  * accumulated, so nothing travels more than a few places from where its
  * priority put it.
  */
  static constexpr std::size_t max_postponed = 8;

  static void run_ordered(std::vector<Element_type>& candidates,
                          Operation& op, C3t3& c3t3)
  {
    tbb::concurrent_queue<Element_type> queue(candidates.begin(), candidates.end());
    tbb::parallel_for(0, tbb::this_task_arena::max_concurrency(),
                      [&](int)
                      {
                        Element_type element;
                        std::vector<Element_type> postponed;
                        while (queue.try_pop(element))
                        {
                          if (!try_apply_one(element, op, c3t3))
                            postponed.push_back(element);

                          if (postponed.size() >= max_postponed)
                            retry_postponed(postponed, op, c3t3);
                        }
                        // Whatever is still held back is taken with the
                        // waiting form: the pass is over for this worker, so
                        // there is nothing else for it to do meanwhile.
                        for (const Element_type& e : postponed)
                          apply_one(e, op, c3t3);
                      });
  }

  // One attempt at each held-back element, in the order they were held back;
  // what still fails stays held back.
  static void retry_postponed(std::vector<Element_type>& postponed,
                              Operation& op, C3t3& c3t3)
  {
    std::size_t kept = 0;
    for (std::size_t i = 0; i < postponed.size(); ++i)
    {
      if (!try_apply_one(postponed[i], op, c3t3))
        postponed[kept++] = postponed[i];
    }
    postponed.resize(kept);
  }

  /**
  * One attempt at `element`. Unlike `apply_one()` it does not wait: a zone it
  * cannot take is given back and the element reported undone, for the caller
  * to come back to.
  */
  static bool try_apply_one(const Element_type& element, Operation& op, C3t3& c3t3)
  {
#ifdef CGAL_TR_LOCKCOUNT
    Lockcount_counters& lc = lockcount_counters();
    ++lc.ops;
#endif
    if(!op.lock_zone(element, c3t3))
    {
#ifdef CGAL_TR_LOCKCOUNT
      ++lc.retries;
#endif
      c3t3.triangulation().unlock_all_elements();
      return false;
    }
    op.execute_operation(element, c3t3);
    c3t3.triangulation().unlock_all_elements();
    return true;
  }

  static std::vector<Element_type>
  gather(tbb::enumerable_thread_specific<std::vector<Element_type>>& per_thread)
  {
    std::vector<Element_type> all;
    std::size_t n = 0;
    for(const auto& v : per_thread) n += v.size();
    all.reserve(n);
    for(const auto& v : per_thread) all.insert(all.end(), v.begin(), v.end());
    return all;
  }

  /**
  * Replays the elements a pass could not take, in parallel rounds.
  *
  * The deferred elements are exactly the ones that CONFLICTED, so they tend to
  * conflict with each other, and replaying them in parallel can replay the
  * same fight: an unconditional round-based replay was measured doing four
  * times the operations on one configuration. The round therefore has to earn
  * the next one -- as soon as a round hands back more than half of what it was
  * given, the rest is done serially, where a zone cannot fail for want of
  * another thread and every element is taken exactly once.
  */
  static void run_deferred(std::vector<Element_type> todo,
                           Operation& op, C3t3& c3t3)
  {
    while(!todo.empty())
    {
      tbb::enumerable_thread_specific<std::vector<Element_type>> again;
      tbb::parallel_for_each(todo,
                             [&](const Element_type& element)
                             {
                               if(!try_apply_one(element, op, c3t3))
                                 again.local().push_back(element);
                             });
      std::vector<Element_type> next = gather(again);
      if(next.empty())
        return;
      if(2 * next.size() > todo.size())   // the round cleared less than half
      {
        for(const Element_type& element : next)
          apply_one(element, op, c3t3);
        return;
      }
      todo.swap(next);
    }
  }

  static void run_unordered(std::vector<Element_type>& candidates,
                            Operation& op, C3t3& c3t3)
  {
    // No shuffle. It was introduced to spread threads over the mesh, and it
    // does not pay for itself: removing it is -2.807% wall time on the frozen
    // Tier-A 24 (CLEAR, 20/24 configs faster, instructions flat at +0.066%),
    // with the gain concentrated on the heavy meshes -- 1146193_cdt_1.5
    // -15.40%, 65617_cdt_0.5 -7.59%. Same work, better order: get_elements()
    // already produces candidates in an order the shuffle was destroying.
    //
    // It was also the largest single source of variance in the measurement
    // rig. Reseeding from random_device on every call makes every run remesh a
    // different sequence, which no replication inside one screen can average
    // out. A/A on the same binary, per-config wall sd falls 7.192% -> 1.401%
    // and the instruction null +0.827% -> -0.043%. One thread becomes
    // deterministic, which is what made an operation-level change measurable
    // at all.
    //
    // A worker that cannot take an element's zone does NOT wait for it. It
    // sets the element aside and takes the next one; what is set aside is
    // replayed once the pass has joined. Waiting was measured as the larger
    // half of the parallel path's cost -- 62-74% of the extra instructions a
    // 4-thread run executes over a 1-thread run are kernel instructions, and
    // they are `sched_yield()` in the retry loop, which almost never finds
    // another runnable task to switch to.
    //
    // Only the UNORDERED operations may do this. An ordered one takes its
    // candidate order from get_elements() for a reason, and deferring there
    // changes the trajectory: deferring every operation was measured running
    // 6.6 M operations against 4.55 M at one thread.
    //
    // At one thread nothing is ever deferred -- `lock_zone()` cannot fail when
    // no other thread holds anything -- so the single-threaded result is
    // exactly what it was.
    tbb::enumerable_thread_specific<std::vector<Element_type>> deferred;
    tbb::parallel_for_each(candidates,
                           [&](const Element_type& element)
                           {
                             if(!try_apply_one(element, op, c3t3))
                               deferred.local().push_back(element);
                           });
    run_deferred(gather(deferred), op, c3t3);
  }
};
#endif // CGAL_LINKED_WITH_TBB

/**
* Selects the execution strategy from the triangulation's concurrency tag.
* `Elementary_operation_executor<Op, Tag>` is the executor to instantiate.
*/
template <typename Operation, typename ConcurrencyTag>
struct Elementary_operation_executor_selector
{
  using type = Elementary_operation_execution_sequential<Operation>;
};

#ifdef CGAL_LINKED_WITH_TBB
template <typename Operation>
struct Elementary_operation_executor_selector<Operation, CGAL::Parallel_tag>
{
  using type = Elementary_operation_execution_parallel<Operation>;
};
#endif

template <typename Operation, typename ConcurrencyTag>
using Elementary_operation_executor =
  typename Elementary_operation_executor_selector<Operation, ConcurrencyTag>::type;

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H
