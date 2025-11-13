#include <CGAL/assertions.h>
#include <CGAL/bisect_failures.h>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <vector>

// Test: bisecting a vector to find which elements cause a "failure"
// In this test, we'll consider the vector "bad" if it contains the
// numbers 42, 99, or 100. We track different failure scenarios:
// - 42 and 99 together: throws std::runtime_error
// - 99 alone: throws std::runtime_error (different message)
// - 42 alone: throws CGAL::Failure_exception via CGAL_error_msg
// - 100 alone: returns exit code 100 (non-exception failure)

using Test_data = std::vector<int>;

auto display_elements = [](const Test_data& data) {
  std::stringstream ss;
  ss << "[ ";
  for(const auto& elem : data) {
    ss << elem << " ";
  }
  ss << "]";
  return ss.str();
};

int test(int test_id,
         const Test_data& data,
         bool expect_failure,
         const Test_data& expected_bad_elements,
         const Test_data& unexpected_bad_elements)
{
  Test_data current_data;
  Test_data other_error_data;
  Test_data final_bad_data;
  Test_data other__should_not_be_used;

  std::cout << "Input data: " << display_elements(data) << "\n";

  // Get size function: returns the number of elements
  auto get_size = [](const Test_data& data) -> std::size_t {
    return data.size();
  };

  // Simplify function: removes elements from [start, end)
  auto simplify = [](Test_data& data, std::size_t start, std::size_t end) -> bool {
    if(start >= data.size()) {
      return false;
    }

    const std::size_t actual_end = (std::min)(end, data.size());
    data.erase(data.begin() + start, data.begin() + actual_end);

    return true;
  };

  // Run function: tests if the data is "bad" (contains 42, 99, or 100)
  auto run = [](const Test_data& data) -> int {
    // Check if the "bad" elements are present
    auto it42 = std::find(data.begin(), data.end(), 42);
    auto it99 = std::find(data.begin(), data.end(), 99);
    auto it100 = std::find(data.begin(), data.end(), 100);

    if(it42 != data.end() && it99 != data.end()) {
      throw std::runtime_error("Found bad elements: 42 and 99");
    }
    if(it99 != data.end()) {
      throw std::runtime_error("Found bad element: 99");
    }
    if(it42 != data.end()) {
      CGAL_error_msg("Found bad element: 42");
    }

    if(it100 != data.end()) {
      return 100;
    }

    return EXIT_SUCCESS;
  };

  // Save function: captures saved data for validation
  auto save = [&](const Test_data& data, std::string name) {
    const std::map<std::string, Test_data*> name_to_data = {
      { "bad", &final_bad_data },
      { "final_bad", &final_bad_data },  // Both map to final_bad_data
      { "error", &other_error_data },
      { "current", &current_data }
    };

    // Don't save files during test
    std::cout << "(fake) Saving data of size " << data.size() << " to file " << name << "\n  "
              << display_elements(data) << "\n";

    auto it = name_to_data.find(name);
    Test_data* target = (it != name_to_data.end()) ? it->second : &other__should_not_be_used;
    *target = data;

  };

  // Run bisection
  int result_code = EXIT_SUCCESS;
  bool caught_exception = false;
  try {
    result_code = CGAL::bisect_failures(data, get_size, simplify, run, save);

    // If we get here, bisection completed without throwing
    // (either no failure, or failure with non-throwing exit code)
  } catch(const std::exception& e) {
    // Expected: bisect_failures re-runs the minimal failing case which throws std::runtime_error
    std::string msg = e.what();
    assert((msg.find("42") != std::string::npos || msg.find("99") != std::string::npos) &&
           "Exception should mention a bad element");
    caught_exception = true;
  }

  // Validate results
  if(expect_failure != (caught_exception || result_code != EXIT_SUCCESS)) {
    std::cerr << "Test failed! Expected failure: " << expect_failure << ", but got "
              << (caught_exception ? "" : "no ") << "exception and result code is " << result_code << ".\n";
    abort();
  }

  if(expected_bad_elements != final_bad_data ||
     unexpected_bad_elements != other_error_data ||
     (!expect_failure && !final_bad_data.empty()) ||
     other__should_not_be_used.size() != 0)
  {
    std::cerr << "Test " << test_id << " failed!";
    std::cerr << "\n  Expected bad elements: ";
    std::cerr << display_elements(expected_bad_elements);
    std::cerr << "\n  Unexpected bad elements: ";
    std::cerr << display_elements(unexpected_bad_elements);
    std::cerr << "\n  Final bad elements: ";
    std::cerr << display_elements(final_bad_data);
    std::cerr << "\n  Other error elements: ";
    std::cerr << display_elements(other_error_data);
    std::cerr << "\n  Other elements (should not be used): ";
    std::cerr << display_elements(other__should_not_be_used);
    abort();
  }

  return EXIT_SUCCESS;
}

int main() {
#if defined(CGAL_NDEBUG)
  std::cerr << "Error: This test requires CGAL assertions to be enabled.\n"
               "       Please compile without NDEBUG and CGAL_NDEBUG.\n";
  return 0;
#endif
  std::cout << "=== Edge Cases ===\n\n";

  std::cout << "## Test 1: Empty data - should succeed\n";
  Test_data data_empty;
  test(1, data_empty, false, {}, {});

  std::cout << "## Test 2: No bad elements - should succeed\n";
  Test_data data_good = {1, 2, 3, 4, 5};
  test(2, data_good, false, {}, {});

  std::cout << "\n=== Single Element Failures ===\n\n";

  std::cout << "## Test 3: Only element 42 (CGAL exception)\n";
  // Element 42 alone triggers CGAL_error_msg
  test(3, {42, 50}, true, {42}, {});

  std::cout << "## Test 4: Only element 99 (std::runtime_error)\n";
  // Element 99 alone triggers std::runtime_error
  test(4, {50, 99}, true, {99}, {});

  std::cout << "## Test 5: Only element 100 (exit code)\n";
  // Element 100 alone causes exit code 100
  test(5, {50, 100}, true, {100}, {});

  std::cout << "\n=== Pair Combinations ===\n\n";

  std::cout << "## Test 6: Elements {42, 99} starting from size 2 - cannot reduce\n";
  // With only 2 elements, bisection can't keep both together
  // Each element is detected separately as different error types
  test(6, {42, 99}, true, {}, {42});

  std::cout << "## Test 7: Elements {42, 99} starting from size 3 - minimal case\n";
  // Starting with 3 elements, algorithm successfully bisects down to {42, 99}
  test(7, {42, 50, 99}, true, {42, 99}, {42});

  std::cout << "## Test 8: Elements {42, 100} (CGAL exception + exit code)\n";
  // Element 42 (CGAL exception) is the main failure
  // Element 100 (exit code) encountered as different error during bisection
  test(8, {42, 100}, true, {42}, {100});

  std::cout << "## Test 9: Elements {99, 100} (std::runtime_error + exit code)\n";
  // Element 99 (std::runtime_error) is the main failure
  // Element 100 (exit code) encountered as different error during bisection
  test(9, {99, 100}, true, {99}, {100});

  std::cout << "\n=== Triple Combination ===\n\n";

  std::cout << "## Test 10: All three elements {42, 99, 100}\n";
  // Algorithm finds {42, 99} as main failure (std::runtime_error)
  // Encounters element 42 alone (CGAL exception) as different error during bisection
  test(10, {42, 99, 100}, true, {42, 99}, {42});

  std::cout << "\n=== Larger Dataset Tests ===\n\n";

  std::cout << "## Test 11: Elements {42, 99} at start of range\n";
  Test_data data_start = {42, 43, 44, 99, 100, 101};
  // Main failure is {42, 99}, encounters element 42 alone during bisection
  test(11, data_start, true, {42, 99}, {42});

  std::cout << "## Test 12: Elements {42, 99} in middle of small dataset\n";
  Test_data data_small = {40, 41, 42, 43, 98, 99, 100, 101};
  // Main failure is {42, 99}, encounters element 42 alone during bisection
  test(12, data_small, true, {42, 99}, {42});

  std::cout << "## Test 13: Large dataset (100 elements) with {42, 99}\n";
  Test_data data_100(100);
  std::iota(data_100.begin(), data_100.end(), 0);
  // Successfully bisects from 100 elements down to {42, 99}
  // Encounters element 42 alone (CGAL exception) during bisection
  test(13, data_100, true, {42, 99}, {42});

  std::cout << "## Test 14: Large sparse dataset (200 elements) with {99}\n";
  Test_data data_200(200);
  std::iota(data_200.begin(), data_200.end(), 50);
  // Successfully bisects from 200 elements down to {99}
  test(14, data_200, true, {99}, {100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149});

  std::cout << "\n✓ All tests passed!\n";
  return 0;
}
