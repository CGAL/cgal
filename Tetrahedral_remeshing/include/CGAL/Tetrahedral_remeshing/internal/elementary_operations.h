#ifndef CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H
#define CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Iterator_range.h>

#include <vector>
#include <string>
#include <atomic>
#include <utility>
#include <iostream>
#include <algorithm>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/task_group.h>
#include <tbb/combinable.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#endif

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

  // Pure element logic methods
  virtual bool should_process_element(const ElementType& e, const C3t3& c3t3) const = 0;
  virtual ElementSource get_element_source(const C3t3& c3t3) const = 0;
  virtual bool can_apply_operation(const ElementType& e, const C3t3& c3t3) const = 0;
  virtual bool lock_zone(const ElementType& e, const C3t3& c3t3) const = 0;
  virtual bool execute_operation(const ElementType& e, C3t3& c3t3) = 0;
  virtual std::string operation_name() const = 0;

  // Pre/post operation hooks
  virtual bool execute_pre_operation(const ElementType& e, C3t3& c3t3) {
    return true; // Default implementation does nothing
  }

  virtual bool execute_post_operation(const ElementType& e, C3t3& c3t3) {
    return true; // Default implementation does nothing
  }
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
  void ensure_lock_data_structure_initialized(C3t3& c3t3) {
#ifdef CGAL_CONCURRENT_TETRAHEDRAL_REMESHING
    auto& triangulation = c3t3.triangulation();
    if (!triangulation.get_lock_data_structure()) {
      // Create a lock data structure using the C3T3's bounding box
      static typename C3t3::Triangulation::Lock_data_structure lock_ds(
        c3t3.bbox(),  // Use C3T3's bounding box
        8  // Default grid size (same as Mesh_3 default)
      );
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
    std::vector<ElementType> candidates;

    // Get all elements from the operation
    auto elements = op.get_element_source(c3t3);

    for (const auto& element : elements) {
      if (op.should_process_element(element, c3t3)) {
        candidates.push_back(element);
      }
    }
    return candidates;
  }

  bool apply_operation_on_elements(
    std::vector<ElementType>& elements,
    Operation& op,
    C3t3& c3t3) override {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Executing operation sequentially: " << op.operation_name() << std::endl;
    #endif

    bool success = false;

    for (const auto& element : elements) {
      if (!op.can_apply_operation(element, c3t3))
        continue;

      if (op.lock_zone(element,c3t3)) {
        // Execute pre-operation phase for this element
        op.execute_pre_operation(element, c3t3);

        if (op.execute_operation(element, c3t3)) {
          success = true;
        }

        // Execute post-operation phase for this element
        op.execute_post_operation(element, c3t3);

      }
      c3t3.triangulation().unlock_all_elements();
    }

    return success;
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

  std::vector<ElementType> collect_candidates(const Operation& op, const C3t3& c3t3) const override {
    // Get all elements from the operation
    auto elements = op.get_element_source(c3t3);

    // Use blocked_range for parallel processing with Iterator_range
    tbb::blocked_range<std::size_t> range(0, elements.size());
    
    // First pass: count elements that pass the filter
    std::atomic<std::size_t> total_count(0);
    tbb::parallel_for(range, [&](const tbb::blocked_range<std::size_t>& r) {
      auto it = elements.begin();
      std::advance(it, r.begin());
      for (std::size_t i = 0; i < r.size(); ++i, ++it) {
        if (op.should_process_element(*it, c3t3)) {
          total_count.fetch_add(1);
        }
      }
    });

    // Pre-allocate result vector
    std::vector<ElementType> candidates;
    candidates.reserve(total_count);

    // Second pass: collect elements using thread-local vectors
    tbb::combinable<std::vector<ElementType>> local_candidates;
    
    tbb::parallel_for(range, [&](const tbb::blocked_range<std::size_t>& r) {
      auto it = elements.begin();
      std::advance(it, r.begin());
      for (std::size_t i = 0; i < r.size(); ++i, ++it) {
        if (op.should_process_element(*it, c3t3)) {
          local_candidates.local().push_back(*it);
        }
      }
    });

    // Combine all local vectors into final result
    local_candidates.combine_each([&](const std::vector<ElementType>& local) {
      candidates.insert(candidates.end(), local.begin(), local.end());
    });

    return candidates;
  }

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

    std::atomic<bool> success(false);

    // TODO: GRAIN_SIZE should be fine tuned
    const size_t GRAIN_SIZE = std::max(elements.size() / (std::thread::hardware_concurrency() * 4), size_t(8));

    tbb::parallel_for(
      tbb::blocked_range<size_t>(0, elements.size(), GRAIN_SIZE),
      [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          const auto& element = elements[i];
          bool affected_zone_locked = op.lock_zone(element, c3t3);
          if (affected_zone_locked) {
            if (!op.can_apply_operation(element, c3t3)) {
              c3t3.triangulation().unlock_all_elements();
              continue;
            }
            if (!op.execute_pre_operation(element, c3t3)) {
              c3t3.triangulation().unlock_all_elements();
              continue;
            }
            if (!op.execute_operation(element, c3t3)) {
              c3t3.triangulation().unlock_all_elements();
              continue;
            }
            if (!op.execute_post_operation(element, c3t3)) {
              c3t3.triangulation().unlock_all_elements();
              continue;
            }
            success.store(true, std::memory_order_relaxed);
          }
          c3t3.triangulation().unlock_all_elements();
        }
      }
    );
    return success.load();
  }
};
} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H