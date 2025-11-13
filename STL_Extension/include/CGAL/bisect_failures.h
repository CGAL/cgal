// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau
//
// The documentation has been partially generated using Github Copilot.

#ifndef CGAL_BISECT_FAILURES_H
#define CGAL_BISECT_FAILURES_H

#include <CGAL/config.h>

#include <CGAL/exceptions.h>
#include <CGAL/utility.h>

#include <cmath>
#include <exception>
#include <iostream>
#include <optional>
#include <string>

namespace CGAL {

/**
 * \ingroup PkgSTLExtensionUtilities
 *
 * \brief Bisects input data by iteratively simplifying it to identify the minimal failing case.
 *
 * This debugging utility helps identify minimal test cases when complex input data causes failures.
 * It works by iteratively simplifying the data and testing whether the failure persists,
 * using a bisection-like approach to narrow down to the smallest failing case.
 *
 * The algorithm divides the input data into "buckets" and systematically removes each bucket
 * to test if the failure persists. It starts with a coarse granularity (ratio=0.5, removing
 * half the elements) and automatically becomes more fine-grained (dividing ratio by 2) when
 * no fault is found. When a failure is found, it restarts the bisection with the smaller
 * failing data, progressively narrowing down to the minimal case.
 *
 * \tparam InputData The type of input data to bisect (must be copyable and assignable)
 * \tparam GetSizeFn Function object type: `std::size_t GetSizeFn(const InputData& data)`
 * \tparam SimplifyFn Function object type: `bool SimplifyFn(InputData& data, std::size_t start, std::size_t end)`
 * \tparam RunFn Function object type: `int RunFn(const InputData& data)`
 * \tparam SaveFn Function object type: `void SaveFn(const InputData& data, const std::string& filename_prefix)`
 *
 * \param data The input data to bisect
 * \param get_size_fn Function that returns the "size" of the data (e.g., number of elements).
 * \param simplify_fn Function that simplifies the data by removing elements from `[start, end)`.
 *                    Should return `true` if simplification succeeded, `false` otherwise.
 * \param run_fn Function that tests the data. Should return 0 (`EXIT_SUCCESS`) on success, non-zero on failure.
 *                May also throw exceptions to indicate failure.
 * \param save_fn Function that saves the data to a file or output. Its second parameter
 *                (`filename_prefix`) indicates the context (e.g., "bad", "final_bad", "error", "current")
 *                and can be used to name the output accordingly.
 *
 * \return Exit code: 0 (EXIT_SUCCESS) if no failures found, non-zero otherwise
 *
 * The algorithm:
 * 1. Tests the full data first to verify it fails and capture the failure pattern
 * 2. Starts with a ratio of 0.5 (removing 50% of elements) and divides data into "buckets"
 * 3. For each bucket, creates a simplified version by removing that portion
 * 4. Tests the simplified version with `run_fn`
 * 5. If it fails with the same pattern as the original, saves it as "bad" and restarts bisection with this smaller dataset
 * 6. If it succeeds or fails differently, tries the next bucket
 * 7. After a complete pass with no matching failures found, reduces the ratio by half (0.5 → 0.25 → 0.125...)
 * 8. Repeats until no further simplification is possible (minimal failing case found)
 * 9. Saves the minimal failing case as "final_bad" and returns its exit code
 *
 * \snippet STL_Extension/bisect_failures.cpp bisect_failures_snippet
 */
template<typename InputData, typename GetSizeFn, typename SimplifyFn, typename RunFn, typename SaveFn>
int bisect_failures(const InputData& data,
                    GetSizeFn get_size_fn,
                    SimplifyFn simplify_fn,
                    RunFn run_fn,
                    SaveFn save_fn)
{
  // Redirect temporarily cout to cerr, and clog to stdout, for debug output
  auto* old_clog_buf = std::clog.rdbuf();
  auto* old_cout_buf = std::cout.rdbuf();
  auto _ = make_scope_exit([&]() {
          std::clog.rdbuf(old_clog_buf);
          std::cout.rdbuf(old_cout_buf);
  });
  std::clog.rdbuf(old_cout_buf);
  std::cout.rdbuf(std::cerr.rdbuf());
  // The following code uses std::clog to display its messages, and the user's
  // code (run_fn) usages of std::cout will go to std::cerr.

  auto do_run = [&](InputData current_data) {
    std::optional<Failure_exception> cgal_exc;
    std::optional<std::string> std_exc_msg;
    int exit_code = EXIT_SUCCESS;
    try {
      exit_code = run_fn(current_data);
    } catch(Failure_exception& e) {
      cgal_exc = e;
      std::clog << "    CAUGHT CGAL EXCEPTION: " << e.what() << '\n';
    } catch(std::exception& e) {
      std::clog << "    CAUGHT EXCEPTION: " << e.what() << '\n';
      std_exc_msg = e.what();
    }
    if(exit_code != EXIT_SUCCESS)
      std::clog << "    RUN RETURNED EXIT CODE: " << exit_code << '\n';
    return std::make_tuple(cgal_exc, std_exc_msg, exit_code);
  };

  std::clog << "First run of the algorithm on full data, to check how it fails\n";
  auto [initial_cgal_exception, initial_std_exception, initial_exit_code] = do_run(data);
  if(!initial_cgal_exception && !initial_std_exception && initial_exit_code == EXIT_SUCCESS) {
    std::clog << "Initial run succeeded, no failure to bisect\n";
    return EXIT_SUCCESS;
  }

  double ratio = 0.5;  // Start with removing half the elements

  InputData bad_data{data};
  InputData working_data{data};

  int exit_code = EXIT_SUCCESS;

  while(true) {
    std::size_t nb_buckets = static_cast<std::size_t>(std::floor(1.0 / ratio)) + 1;
    std::clog << "RATIO: " << ratio << '\n';
    bool found_fault_this_pass = false;

    std::size_t nb_to_skip = 0;
    for(std::size_t bucket = 0; bucket < nb_buckets;) {
      const auto data_size = get_size_fn(working_data);
      nb_to_skip = static_cast<std::size_t>(std::round(data_size * ratio));
      if(nb_to_skip < 1) {
        nb_to_skip = 1;
        nb_buckets = data_size;
      }

      const auto start = (std::min)(bucket * nb_to_skip, data_size);
      const auto end = (std::min)(start + nb_to_skip, data_size);
      std::clog << "  SKIP from " << start << " to " << end << '\n';


      // Try to simplify the data
      if(simplify_fn(working_data, start, end)) {
        const auto new_size = get_size_fn(working_data);
        std::clog << "    size after simplification: " << new_size << '\n';

        if(new_size >= data_size) {
          std::clog << "    ERROR: could not simplify data\n";
          working_data = bad_data;
          ++bucket;
          continue;
        }

        // Save current state
        save_fn(working_data, "current");

        auto [cgal_exception, std_exception, this_run_exit_code] = do_run(working_data);

        bool same_exception = false;
        bool same_exit_code = false;

        if(cgal_exception) {
          if(initial_cgal_exception &&
             cgal_exception->message() == initial_cgal_exception->message() &&
             cgal_exception->expression() == initial_cgal_exception->expression() &&
             cgal_exception->library() == initial_cgal_exception->library() &&
             cgal_exception->filename() == initial_cgal_exception->filename() &&
             cgal_exception->line_number() == initial_cgal_exception->line_number())
          {
            same_exception = true;
          }
        } else if(std_exception) {
          // Check if this is the same type of failure we're looking for
          if(initial_std_exception &&
             *initial_std_exception == *std_exception)
          {
            same_exception = true;
          }
        } else if(this_run_exit_code != EXIT_SUCCESS) {
          if(this_run_exit_code == initial_exit_code) {
            same_exit_code = true;
          }
        }
        if(this_run_exit_code != EXIT_SUCCESS) {
          exit_code = this_run_exit_code;
        }
        if(same_exception || same_exit_code) {
            std::clog << "    -> BAD DATA! (size: " << get_size_fn(working_data) << ")\n";
            save_fn(working_data, "bad");
            bad_data = working_data;
            found_fault_this_pass = true;
            bucket = 0; // Reset to bisect further
            continue;
        }
        if (cgal_exception || std_exception || this_run_exit_code != EXIT_SUCCESS) {
            // Different type of error - log it but continue
            if(exit_code == EXIT_SUCCESS) exit_code = EXIT_FAILURE;
            std::clog << "    -> ERROR DATA (different error type)\n";
            save_fn(working_data, "error");
            std::clog << "       go on...\n";
          } else {
            std::clog << "    -> GOOD DATA :-( (size: " << get_size_fn(working_data) << ")\n";
          }
        }

      // Reset to bad_data for next iteration
      working_data = bad_data;
      ++bucket;
    }

    // After completing a full pass through all buckets
    if(!found_fault_this_pass) {
      if(nb_to_skip <= 1) {
        // Cannot subdivide further - we've found the minimal failing case
        break;
      }
      // No fault found at this ratio - make ratio smaller (more granular)
      ratio = ratio / 2.0;

      nb_buckets = static_cast<std::size_t>(std::floor(1.0 / ratio)) + 1;
      std::clog << "  No fault found at this ratio. Reducing ratio to: " << ratio << '\n';
    }
  }

  if(get_size_fn(bad_data) < get_size_fn(data)) {
    std::clog << "FINAL BAD DATA: " << get_size_fn(bad_data) << " elements\n";
    save_fn(bad_data, "final_bad");
    return run_fn(bad_data);
  }

  return exit_code;
}

} // namespace CGAL

#endif // CGAL_BISECT_FAILURES_H
