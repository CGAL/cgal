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

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/task_group.h>
#include <tbb/combinable.h>
#endif

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

// Generic lock manager that works with any lockable element type
template<typename LockableHandle>
class LockManager {
public:
  using Lock_zone = std::vector<LockableHandle>;

  LockManager() = default;

  bool try_lock_zone(const Lock_zone& zone) {
    // First try to acquire all locks
    //for (const LockableHandle& h : zone) {
    //  if (!h->try_lock()) {
    //    // If any lock fails, release all acquired locks and return false
    //    unlock_zone(zone);
    //    return false;
    //  }
    //}
    //return true;
    return false;
  }

  void unlock_zone(const Lock_zone& zone) {
    //for (const LockableHandle& h : zone) {
    //  h->unlock();
    //}
  }
};

// Base class for elementary operations
template<typename C3t3_,typename ElementIteratorType_,typename LockElementType_>
class ElementaryOperation {
public:
  using C3t3 = C3t3_;
  using Triangulation = typename C3t3::Triangulation;
  using ElementIteratorType = ElementIteratorType_;
  using LockElementType = LockElementType_;
  using Lock_zone = std::vector<LockElementType>;

  ElementaryOperation() = default;
  virtual ~ElementaryOperation() = default;

  // Returns true if the element should be processed
  virtual bool should_process_element(const ElementIteratorType& e, const C3t3& c3t3) const = 0;

  // Returns the range of elements to process
  virtual typename CGAL::Iterator_range<ElementIteratorType>
  get_element_iterators(const C3t3& c3t3) const = 0;

  // Returns true if the operation can be applied to the element
  virtual bool can_apply_operation(const ElementIteratorType& e, const C3t3& c3t3) const = 0;

  // Returns the zone that needs to be locked for this operation
  virtual Lock_zone get_lock_zone(const ElementIteratorType& e, const C3t3& c3t3) const = 0;

  // Pre-operation setup, returns false if operation should be aborted
  virtual bool execute_pre_operation(const ElementIteratorType& e, C3t3& c3t3) = 0;

  // Execute the operation, returns false if operation failed
  virtual bool execute_operation(const ElementIteratorType& e, C3t3& c3t3) = 0;

  // Post-operation cleanup, returns false if operation should be rolled back
  virtual bool execute_post_operation(const ElementIteratorType& e, C3t3& c3t3) = 0;

  // Returns the name of the operation for debugging/logging
  virtual std::string operation_name() const = 0;
};

// Base class for operation execution strategies
template<typename Operation>
class ElementaryOperationExecution {
public:
  using C3t3 = typename Operation::C3t3;
  using Cell_handle = typename C3t3::Triangulation::Cell_handle;
  using ElementIteratorType = typename Operation::ElementIteratorType;
  ElementaryOperationExecution() = default;
  ElementaryOperationExecution(const ElementaryOperationExecution&) = default;
  ElementaryOperationExecution(ElementaryOperationExecution&&) = default;
  ElementaryOperationExecution& operator=(const ElementaryOperationExecution&) = default;
  ElementaryOperationExecution& operator=(ElementaryOperationExecution&&) = default;
  virtual ~ElementaryOperationExecution() = default;

  virtual std::vector<ElementIteratorType> collect_candidates(
    const Operation& op, const C3t3& c3t3) const = 0;

  virtual bool apply_operation_on_elements(
    std::vector<ElementIteratorType>& elements,
    Operation& op,
    C3t3& c3t3,
    LockManager<Cell_handle>& lock_manager) = 0;

  void execute(Operation& op, C3t3& c3t3) {
    LockManager<Cell_handle> lock_manager;
    auto candidates = collect_candidates(op, c3t3);
    apply_operation_on_elements(candidates, op, c3t3, lock_manager);
  }
};

// Sequential execution implementation
template<typename Operation>
class ElementaryOperationExecutionSequential : public ElementaryOperationExecution<Operation> {
public:
  using Base = ElementaryOperationExecution<Operation>;
  using C3t3 = typename Operation::C3t3;
  using ElementIteratorType = typename Operation::ElementIteratorType;
  using Cell_handle = typename C3t3::Triangulation::Cell_handle;

  ElementaryOperationExecutionSequential() = default;
  ElementaryOperationExecutionSequential(const ElementaryOperationExecutionSequential&) = default;
  ElementaryOperationExecutionSequential(ElementaryOperationExecutionSequential&&) = default;
  ElementaryOperationExecutionSequential& operator=(const ElementaryOperationExecutionSequential&) = default;
  ElementaryOperationExecutionSequential& operator=(ElementaryOperationExecutionSequential&&) = default;

  std::vector<ElementIteratorType> collect_candidates(
    const Operation& op, const C3t3& c3t3) const override 
  {
    std::vector<ElementIteratorType> candidates;
    auto [begin, end] = op.get_element_iterators(c3t3);
    for (auto it = begin; it != end; ++it) {
      if (op.should_process_element(it, c3t3)) {
        candidates.push_back(it);
      }
    }
    return candidates;
  }

  bool apply_operation_on_elements(
    std::vector<ElementIteratorType>& elements,
    Operation& op,
    C3t3& c3t3,
    LockManager<Cell_handle>& lock_manager) override 
  {
    bool success = false;
    for (const ElementIteratorType& e : elements) {
      if (!op.can_apply_operation(e, c3t3))
        continue;

      auto lock_zone = op.get_lock_zone(e, c3t3);
      if (lock_manager.try_lock_zone(lock_zone)) {
        if (op.execute_pre_operation(e, c3t3)) {
          if (op.execute_operation(e, c3t3)) {
            op.execute_post_operation(e, c3t3);
            success = true;
          }
        }
        lock_manager.unlock_zone(lock_zone);
      }
    }
    return success;
  }

};

#ifdef CGAL_LINKED_WITH_TBB
// Parallel execution implementation
template<typename Operation>
class ElementaryOperationExecutionParallel : public ElementaryOperationExecution<Operation> {
public:
  using Base = ElementaryOperationExecution<Operation>;
  using C3t3 = typename Operation::C3t3;
  using ElementIteratorType = typename Operation::ElementIteratorType;
  using Cell_handle = typename C3t3::Triangulation::Cell_handle;

  ElementaryOperationExecutionParallel() = default;
  ElementaryOperationExecutionParallel(const ElementaryOperationExecutionParallel&) = default;
  ElementaryOperationExecutionParallel(ElementaryOperationExecutionParallel&&) = default;
  ElementaryOperationExecutionParallel& operator=(const ElementaryOperationExecutionParallel&) = default;
  ElementaryOperationExecutionParallel& operator=(ElementaryOperationExecutionParallel&&) = default;

  std::vector<ElementIteratorType> collect_candidates(
    const Operation& op, const C3t3& c3t3) const override 
  {
    tbb::combinable<std::vector<ElementIteratorType>> local_candidates;
    auto [begin, end] = op.get_element_iterators(c3t3);
    
    tbb::parallel_for(tbb::blocked_range<ElementIteratorType*>(begin, end),
      [&](const tbb::blocked_range<ElementIteratorType*>& r) {
        auto& my_candidates = local_candidates.local();
        for (auto it = r.begin(); it != r.end(); ++it) {
          if (op.should_process_element(*it, c3t3)) {
            my_candidates.push_back(*it);
          }
        }
      });
    
    std::vector<ElementIteratorType> candidates;
    local_candidates.combine_each([&](const std::vector<ElementIteratorType>& local) {
      candidates.insert(candidates.end(), local.begin(), local.end());
    });
    return candidates;
  }

  bool apply_operation_on_elements(
    std::vector<ElementIteratorType>& elements,
    Operation& op,
    C3t3& c3t3,
    LockManager<Cell_handle>& lock_manager) override 
  {
    tbb::task_group group;
    std::atomic<bool> success{false};

    for (const ElementIteratorType& e : elements) {
      group.run([&]() {
        if (!op.can_apply_operation(e, c3t3))
          return;

        auto lock_zone = op.get_lock_zone(e, c3t3);
        if (lock_manager.try_lock_zone(lock_zone)) {
          if (op.execute_pre_operation(e, c3t3)) {
            if (op.execute_operation(e, c3t3)) {
              op.execute_post_operation(e, c3t3);
              success = true;
            }
          }
          lock_manager.unlock_zone(lock_zone);
        }
      });
    }

    group.wait();
    return success;
  }
};
#endif // CGAL_LINKED_WITH_TBB

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_OPERATIONS_H 