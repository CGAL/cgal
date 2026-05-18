#include <CGAL/IO/Indenting_ostream.h>
#include <iostream>
#include <sstream>

/*!
 * \example Stream_support/indenting_ostream.cpp
 *
 * This example demonstrates how to use CGAL::IO::Basic_indenting_streambuf
 * and CGAL::IO::basic_indenting_stream_guard for formatting output with
 * automatic indentation.
 */

void demonstrate_basic_usage() {
  std::cout << "=== Basic Indenting Usage ===\n\n";

  // Example 1: Basic indentation with std::cout
  {
    std::cout << "1. Basic indentation with std::cout:\n";
    CGAL::IO::Indenting_stream_guard guard(std::cout, "  ");

    std::cout << "This line is indented\n";
    std::cout << "So is this one\n";
    std::cout << "Multi-line output:\n";
    std::cout << "Line 1\nLine 2\nLine 3\n";
  }
  std::cout << "Back to normal indentation\n\n";

  // Example 2: Nested indentation levels
  {
    std::cout << "2. Nested indentation levels:\n";
    CGAL::IO::Indenting_stream_guard level1(std::cout, "| ");

    std::cout << "Level 1 indentation\n";
    {
      CGAL::IO::Indenting_stream_guard level2(std::cout, "| ");
      std::cout << "Level 2 indentation\n";
      {
        CGAL::IO::Indenting_stream_guard level3(std::cout, "| ");
        std::cout << "Level 3 indentation\n";
      }
      std::cout << "Back to level 2\n";
    }
    std::cout << "Back to level 1\n";
  }
  std::cout << "Back to normal\n\n";

  // Example 3: Custom indentation string
  {
    std::cout << "3. Custom indentation string:\n";
    CGAL::IO::Indenting_stream_guard guard(std::cout, ">>> ");

    std::cout << "Custom prefix on each line\n";
    std::cout << "Works with multiple lines too\n";
  }
  std::cout << "Normal output again\n\n";
}

void demonstrate_stringstream_usage() {
  std::cout << "=== Stringstream Usage ===\n\n";

  // Using with stringstream for formatted output
  // Note: We create a custom ostream using the indenting streambuf
  std::ostringstream backing_stream;
  CGAL::IO::Indenting_streambuf indent_buf(*backing_stream.rdbuf(), "    ");
  std::ostream oss(&indent_buf);

  oss << "This goes to stringstream\n";
  oss << "With indentation\n";
  oss << "Multiple lines\n";

  std::cout << "Stringstream contents:\n" << backing_stream.str() << "\n";
}

void simulate_debug_output() {
  std::cout << "=== Simulated Debug Output ===\n\n";

  std::cout << "Algorithm: Starting geometric computation\n";

  {
    CGAL::IO::Indenting_stream_guard guard(std::cout, "  ");
    std::cout << "Phase 1: Input validation\n";
    std::cout << "- Checking vertices: OK\n";
    std::cout << "- Checking faces: OK\n";

    {
      CGAL::IO::Indenting_stream_guard inner_guard(std::cout, "    ");
      std::cout << "Detailed validation:\n";
      std::cout << "- Manifold check: PASSED\n";
      std::cout << "- Orientation check: PASSED\n";
      std::cout << "- Degeneracy check: PASSED\n";
    }

    std::cout << "Phase 1 complete\n\n";

    std::cout << "Phase 2: Triangulation\n";
    std::cout << "- Inserting vertices...\n";
    std::cout << "- Computing Delaunay triangulation...\n";
    std::cout << "- Applying constraints...\n";

    {
      CGAL::IO::Indenting_stream_guard constraint_guard(std::cout, "    ");
      std::cout << "Constraint processing:\n";
      for(int i = 0; i < 3; ++i) {
        std::cout << "- Processing constraint " << i << ": OK\n";
      }
    }

    std::cout << "Phase 2 complete\n";
  }

  std::cout << "Algorithm: Geometric computation finished successfully\n\n";
}

void demonstrate_multiple_streams() {
  std::cout << "=== Multiple Streams with make_indenting_guards ===\n\n";

  std::cout << "Before indentation:\n";
  std::cout << "  cout: Normal output\n";
  std::cerr << "  cerr: Normal output\n";

  // Apply indentation to both streams simultaneously
  {
    std::cout << "\nWith make_indenting_guards (2 spaces):\n";
    auto guards = CGAL::IO::make_indenting_guards(2, std::cout, std::cerr);

    std::cout << "This line is indented on cout\n";
    std::cout << "And this one too\n";

    std::cerr << "This line is indented on cerr\n";
    std::cerr << "And this one too\n";

    // Nested indentation works too
    {
      auto nested_guards = CGAL::IO::make_indenting_guards(2, std::cout, std::cerr);
      std::cout << "Double indented on cout\n";
      std::cerr << "Double indented on cerr\n";
    }

    std::cout << "Back to single indent\n";
    std::cerr << "Back to single indent\n";
  }

  std::cout << "\nAfter indentation scope:\n";
  std::cout << "  cout: Normal output again\n";
  std::cerr << "  cerr: Normal output again\n";

  // Using custom indent string with multiple streams
  {
    std::cout << "\nUsing custom indent string \">> \" with std::cout and std::clog:\n";
    auto guards = CGAL::IO::make_indenting_guards(">> ", std::cout, std::clog);

    std::cout << "cout: Output with custom indent\n";
    std::clog << "clog: Log message with custom indent\n";
    std::cout << "cout: More output\n";
  }

  std::cout << "\n";
}

int main() {
  demonstrate_basic_usage();
  demonstrate_stringstream_usage();
  simulate_debug_output();
  demonstrate_multiple_streams();

  return 0;
}