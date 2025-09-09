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
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Iterator_range.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/task_group.h>
#include <tbb/combinable.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for_each.h>
#include <tbb/concurrent_queue.h>
#endif

#include <atomic>
#include <algorithm>
#include <random>
#include <vector>
#include <thread>
#include <string>
#include <utility>
#include <iostream>
#include <fstream>

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

template<typename C3t3_, typename ElementType_, typename ElementSource_, typename LockElementType_>
class ElementaryOperation {
public:
  using C3t3 = C3t3_;
  using Triangulation = typename C3t3::Triangulation;
  using ElementType = ElementType_;
  using ElementSource = ElementSource_;
  using LockElementType = LockElementType_;
  using Lock_zone = std::vector<LockElementType>;

  ElementaryOperation() = default;
  virtual ~ElementaryOperation() = default;

  // Processing strategy trait - determines execution approach
  virtual bool requires_ordered_processing() const { return false; }
  // Pure element logic methods
  virtual ElementSource get_element_source(const C3t3& c3t3) const = 0;
  virtual bool lock_zone(const ElementType& e, const C3t3& c3t3) const = 0;
  virtual bool execute_operation(const ElementType& e, C3t3& c3t3) = 0;
  virtual std::string operation_name() const = 0;
};

// Base class for operation execution strategies
template<typename Operation>
class ElementaryOperationExecution {
public:
  using C3t3 = typename Operation::C3t3;
  using Cell_handle = typename C3t3::Triangulation::Cell_handle;
  using ElementType = typename Operation::ElementType;

  ElementaryOperationExecution() = default;
  ElementaryOperationExecution(const ElementaryOperationExecution&) = default;
  ElementaryOperationExecution(ElementaryOperationExecution&&) = default;
  ElementaryOperationExecution& operator=(const ElementaryOperationExecution&) = default;
  ElementaryOperationExecution& operator=(ElementaryOperationExecution&&) = default;
  virtual ~ElementaryOperationExecution() = default;

  virtual std::vector<ElementType> collect_candidates(const Operation& op, const C3t3& c3t3) const = 0;

  virtual bool apply_operation_on_elements(
    std::vector<ElementType>& elements,
    Operation& op,
    C3t3& c3t3) = 0;

  // Main execution method that orchestrates the operation
  virtual bool execute(Operation& op, C3t3& c3t3) {
    // Collect candidate elements
    auto candidates = collect_candidates(op, c3t3);

    if (candidates.empty()) {
      return false;
    }

    // Ensure lock data structure is initialized for parallel mode
    ensure_lock_data_structure_initialized(c3t3);

    // Apply operations on elements
    return apply_operation_on_elements(candidates, op, c3t3);
  }

private:
  // Ensure the lock data structure is initialized when needed
  void ensure_lock_data_structure_initialized(C3t3& c3t3)
  {
#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
    auto& triangulation = c3t3.triangulation();
    if (!triangulation.get_lock_data_structure()) {
      // Create a lock data structure using the C3T3's bounding box

        #ifndef USE_TAG_NON_BLOCKING
      static typename C3t3::Triangulation::Lock_data_structure lock_ds(
        c3t3.bbox(),  // Use C3T3's bounding box
        8  // Default grid size (same as Mesh_3 default)
      );
      #else
      static CGAL::Spatial_lock_grid_3<Tag_non_blocking> lock_ds(c3t3.bbox(), 8);
          #endif
      triangulation.set_lock_data_structure(&lock_ds);
    }
#else
    // In sequential mode, no lock data structure is needed
    (void)c3t3;  // Suppress unused parameter warning
#endif
  }
};

// Sequential execution implementation
template<typename Operation>
class ElementaryOperationExecutionSequential : public ElementaryOperationExecution<Operation> {
public:
  using Base = ElementaryOperationExecution<Operation>;
  using C3t3 = typename Operation::C3t3;
  using ElementType = typename Operation::ElementType;
  using Cell_handle = typename C3t3::Triangulation::Cell_handle;

public:
  ElementaryOperationExecutionSequential() = default;
  ElementaryOperationExecutionSequential(const ElementaryOperationExecutionSequential&) = default;
  ElementaryOperationExecutionSequential(ElementaryOperationExecutionSequential&&) = default;
  ElementaryOperationExecutionSequential& operator=(const ElementaryOperationExecutionSequential&) = default;
  ElementaryOperationExecutionSequential& operator=(ElementaryOperationExecutionSequential&&) = default;

  std::vector<ElementType> collect_candidates(const Operation& op, const C3t3& c3t3) const override {
    auto elements = op.get_element_source(c3t3);
    std::vector<ElementType> candidates;
    //candidates.reserve(elements.size());
    //std::copy(elements.begin(), elements.end(), std::back_inserter(candidates));

    if constexpr(std::is_same<decltype(elements), std::vector<ElementType>>::value) {
      candidates = std::move(elements); // directly move if vector
    } else {
      candidates.reserve(std::distance(elements.begin(), elements.end()));
      std::copy(elements.begin(), elements.end(), std::back_inserter(candidates));
    }
    return candidates;
  }

  bool apply_operation_on_elements(
    std::vector<ElementType>& elements,
    Operation& op,
    C3t3& c3t3) override {
    size_t num_successful_locks = 0;
    bool success = true;
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Executing operation sequentially: " << op.operation_name() << "...";
    #endif

    for (const auto& element : elements) {
      if(op.execute_operation(element, c3t3)) {
      num_successful_locks++;
      }

    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << " done num_elements:" << elements.size()
              << ", num_successful_locks : " <<num_successful_locks << std::endl;
    #endif

    return true;
  }
};

template<typename Operation>
class ElementaryOperationExecutionParallel : public ElementaryOperationExecution<Operation> {
public:
  using Base = ElementaryOperationExecution<Operation>;
  using C3t3 = typename Operation::C3t3;
  using ElementType = typename Operation::ElementType;
  using Cell_handle = typename C3t3::Triangulation::Cell_handle;

private:
  size_t m_batch_size;

public:
  ElementaryOperationExecutionParallel(size_t batch_size = 64) : m_batch_size(batch_size) {}
  ElementaryOperationExecutionParallel(const ElementaryOperationExecutionParallel&) = default;
  ElementaryOperationExecutionParallel(ElementaryOperationExecutionParallel&&) = default;
  ElementaryOperationExecutionParallel& operator=(const ElementaryOperationExecutionParallel&) = default;
  ElementaryOperationExecutionParallel& operator=(ElementaryOperationExecutionParallel&&) = default;

  std::vector<ElementType> collect_candidates(const Operation& op, const C3t3& c3t3) const override
  {
    auto elements = op.get_element_source(c3t3);
    std::vector<ElementType> candidates;
    if constexpr(std::is_same<decltype(elements), std::vector<ElementType>>::value) {
      candidates = std::move(elements); // directly move if vector
    } else {
#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
      tbb::combinable<std::vector<ElementType>> local_candidates;
      tbb::parallel_for_each(elements.begin(), elements.end(),
                             [&](const ElementType& element) { local_candidates.local().push_back(element); });
      local_candidates.combine_each([&](const std::vector<ElementType>& local) {
        candidates.insert(candidates.end(), local.begin(), local.end());
      });
#endif
    }
    return candidates;
  }

private:
  // Ordered processing using concurrent queue (for operations like EdgeSplit)
  bool apply_ordered_processing(
    std::vector<ElementType>& elements,
    Operation& op,
    C3t3& c3t3)
    {
#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
    // Create concurrent priority queue for all elements
    tbb::concurrent_queue<ElementType> work_queue(elements.begin(), elements.end());


#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::atomic<size_t> num_successful_locks = 0;
    std::atomic<size_t> num_failed_locks = 0;
#endif
    //size_t num_threads = 8; // Start with single thread for ordered processing
    size_t num_threads = std::thread::hardware_concurrency() / 2;

    // Parallel work stealing from priority queue
    tbb::parallel_for(tbb::blocked_range<size_t>(0, num_threads),
      [&](const tbb::blocked_range<size_t>& r) {
        for (size_t thread_id = r.begin(); thread_id != r.end(); ++thread_id) {
          while (!work_queue.empty()) {
            ElementType element;
            if(!work_queue.try_pop(element)) {
              continue;
            }
            // Process element - priority order maintained by queue
            // Retry until lock is acquired
            bool lock_acquired = false;
            while (!lock_acquired) {
              lock_acquired = op.lock_zone(element, c3t3);
              if (!lock_acquired) {
                // Unlock all elements before retrying
                c3t3.triangulation().unlock_all_elements();
                // Optional: small delay to reduce contention
                std::this_thread::yield();
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
                num_failed_locks++;
#endif
              }
            }
            // Execute operation once lock is acquired
            bool success = op.execute_operation(element, c3t3);
            #ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
            num_successful_locks++;
            #endif
            c3t3.triangulation().unlock_all_elements();
          }
        }
      }
    );

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    {
      const char* csv_path = "lock_stats.csv";
      bool need_header = true;
      {
        std::ifstream ifs(csv_path);
        need_header = !ifs.good() || (ifs.peek() == std::ifstream::traits_type::eof());
      }
      std::ofstream ofs(csv_path, std::ios::app);
      if (need_header) {
        ofs << "operation,num_successful_locks,num_failed_locks,lock_success_rate(%),operation_exec_counter" << std::endl;
      }
//      size_t operation_exec_counter = (op.operation_name() == std::string("Edge Split")) ? g_edge_split_exec_counter : 0;
//      const size_t attempts = static_cast<size_t>(num_successful_locks + num_failed_locks);
//      const double lock_success_rate = (attempts == 0) ? 0.0 : 100.0 * static_cast<double>(num_successful_locks) / //static_cast<double>(attempts);
//      ofs << op.operation_name() << "," << num_successful_locks << "," << num_failed_locks << "," << lock_success_rate << /"," /<< operation_exec_counter << std::endl;
    }
    #endif
#endif //concurrent

    return true;
  }

  // Unordered processing using parallel_for_each (for operations like VertexSmooth, EdgeFlip)
  bool apply_unordered_processing(
    std::vector<ElementType>& elements,
    Operation& op,
    C3t3& c3t3)
    {
#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(elements.begin(), elements.end(), g);

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::atomic<size_t> num_successful_locks = 0;
    std::atomic<size_t> num_failed_locks = 0;
#endif
    tbb::parallel_for_each(elements, [&](const ElementType& element) {
      bool lock_acquired = false;
      while(!lock_acquired) {
        lock_acquired = op.lock_zone(element, c3t3);
        if(!lock_acquired) {
          c3t3.triangulation().unlock_all_elements();
          std::this_thread::yield(); // backoff to reduce contention
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
          num_failed_locks++;
#endif
        }
      }

      bool success = op.execute_operation(element, c3t3);
      #ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      num_successful_locks++;
      #endif
      c3t3.triangulation().unlock_all_elements();
    });
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    {
      const char* csv_path = "lock_stats.csv";
      bool need_header = true;
      {
        std::ifstream ifs(csv_path);
        need_header = !ifs.good() || (ifs.peek() == std::ifstream::traits_type::eof());
      }
      std::ofstream ofs(csv_path, std::ios::app);
      if (need_header) {
        ofs << "operation,num_successful_locks,num_failed_locks,lock_success_rate(%),operation_exec_counter" << std::endl;
      }
//      size_t operation_exec_counter = (op.operation_name() == std::string("Edge Split")) ? g_edge_split_exec_counter : 0;
//      const size_t attempts = static_cast<size_t>(num_successful_locks + num_failed_locks);
//      const double lock_success_rate = (attempts == 0) ? 0.0 : 100.0 * static_cast<double>(num_successful_locks) / static_cast<double>(attempts);
//      ofs << op.operation_name() << "," << num_successful_locks << "," << num_failed_locks << "," << lock_success_rate << "," << operation_exec_counter << std::endl;
    }
    #endif
#endif //concurrent
    return true;
  }

public:
  bool apply_operation_on_elements(
    std::vector<ElementType>& elements,
    Operation& op,
    C3t3& c3t3) override {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Executing operation in parallel: " << op.operation_name() << std::endl;
#endif

    if (elements.empty()) {
      return false;
    }

    // Choose processing strategy based on operation requirements
    if (op.requires_ordered_processing()) {
      return apply_ordered_processing(elements, op, c3t3);
    } else {
      return apply_unordered_processing(elements, op, c3t3);
    }
  }
}; // class ElementaryOperationExecution

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H