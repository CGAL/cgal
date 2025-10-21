#include <CGAL/IO/Color_ostream.h>
#include <iostream>

/*!
 * \example Stream_support/color_ostream.cpp
 *
 * This example demonstrates how to use CGAL::Basic_color_streambuf
 * and CGAL::Basic_color_stream_guard for formatting output with
 * ANSI colors.
 */

void demonstrate_basic_colors() {
  std::cout << "=== Basic Color Usage ===\n\n";

  // Example 1: Single color
  {
    std::cout << "1. Single color (Red):\n";
    CGAL::Color_stream_guard guard(std::cout, CGAL::Ansi_color::Red);

    std::cout << "This text is red\n";
    std::cout << "All lines in this scope are red\n";
  }
  std::cout << "Back to normal color\n\n";

  // Example 2: Different colors
  {
    std::cout << "2. Different colors:\n";

    {
      CGAL::Color_stream_guard green(std::cout, CGAL::Ansi_color::Green);
      std::cout << "Green text\n";
    }

    {
      CGAL::Color_stream_guard blue(std::cout, CGAL::Ansi_color::Blue);
      std::cout << "Blue text\n";
    }

    {
      CGAL::Color_stream_guard yellow(std::cout, CGAL::Ansi_color::Yellow);
      std::cout << "Yellow text\n";
    }
  }
  std::cout << "Back to normal\n\n";

  // Example 3: Bright colors
  {
    std::cout << "3. Bright colors:\n";

    {
      CGAL::Color_stream_guard bright_red(std::cout, CGAL::Ansi_color::BrightRed);
      std::cout << "Bright red text\n";
    }

    {
      CGAL::Color_stream_guard bright_green(std::cout, CGAL::Ansi_color::BrightGreen);
      std::cout << "Bright green text\n";
    }

    {
      CGAL::Color_stream_guard bright_cyan(std::cout, CGAL::Ansi_color::BrightCyan);
      std::cout << "Bright cyan text\n";
    }
  }
  std::cout << "Back to normal\n\n";
}

void demonstrate_combined_colors() {
  std::cout << "=== Combined Colors (Bold, Underline, etc.) ===\n\n";

  // Example 1: Bold text
  {
    std::cout << "1. Bold colors:\n";
    CGAL::Color_stream_guard guard(std::cout,
      std::vector<CGAL::Ansi_color>{CGAL::Ansi_color::Bold, CGAL::Ansi_color::Red});

    std::cout << "This text is bold and red\n";
  }
  std::cout << "Back to normal\n\n";

  // Example 2: Underlined text
  {
    std::cout << "2. Underlined colors:\n";
    CGAL::Color_stream_guard guard(std::cout,
      std::vector<CGAL::Ansi_color>{CGAL::Ansi_color::Underline, CGAL::Ansi_color::Blue});

    std::cout << "This text is underlined and blue\n";
  }
  std::cout << "Back to normal\n\n";

  // Example 3: Bold + underline + color
  {
    std::cout << "3. Multiple attributes:\n";
    CGAL::Color_stream_guard guard(std::cout,
      std::vector<CGAL::Ansi_color>{CGAL::Ansi_color::Bold, CGAL::Ansi_color::Underline, CGAL::Ansi_color::Green});

    std::cout << "This text is bold, underlined, and green\n";
  }
  std::cout << "Back to normal\n\n";
}

void demonstrate_background_colors() {
  std::cout << "=== Background Colors ===\n\n";

  {
    std::cout << "1. Red background:\n";
    CGAL::Color_stream_guard guard(std::cout, CGAL::Ansi_color::BgRed);
    std::cout << "Text with red background\n";
  }
  std::cout << "Back to normal\n\n";

  {
    std::cout << "2. Foreground + background:\n";
    CGAL::Color_stream_guard guard(std::cout,
      std::vector<CGAL::Ansi_color>{CGAL::Ansi_color::Yellow, CGAL::Ansi_color::BgBlue});
    std::cout << "Yellow text on blue background\n";
  }
  std::cout << "Back to normal\n\n";

  {
    std::cout << "3. Bold white on bright green background:\n";
    CGAL::Color_stream_guard guard(std::cout,
      std::vector<CGAL::Ansi_color>{CGAL::Ansi_color::Bold, CGAL::Ansi_color::White, CGAL::Ansi_color::BgBrightGreen});
    std::cout << "Bold white on bright green\n";
  }
  std::cout << "Back to normal\n\n";
}

void demonstrate_nested_colors() {
  std::cout << "=== Nested Color Scopes ===\n\n";

  {
    CGAL::Color_stream_guard outer(std::cout, CGAL::Ansi_color::Blue);
    std::cout << "Outer scope: blue text\n";

    {
      CGAL::Color_stream_guard inner(std::cout, CGAL::Ansi_color::Red);
      std::cout << "Inner scope: red text\n";

      {
        CGAL::Color_stream_guard innermost(std::cout, CGAL::Ansi_color::Green);
        std::cout << "Innermost scope: green text\n";
      }

      std::cout << "Back to inner scope: red text\n";
    }

    std::cout << "Back to outer scope: blue text\n";
  }
  std::cout << "Back to normal\n\n";
}

void simulate_colored_debug_output() {
  std::cout << "=== Simulated Colored Debug Output ===\n\n";

  std::cout << "Algorithm: Starting computation\n";

  {
    CGAL::Color_stream_guard info(std::cout, CGAL::Ansi_color::Cyan);
    std::cout << "INFO: Initializing data structures\n";
  }

  {
    CGAL::Color_stream_guard success(std::cout,
      std::vector<CGAL::Ansi_color>{CGAL::Ansi_color::Bold, CGAL::Ansi_color::Green});
    std::cout << "SUCCESS: Data loaded successfully\n";
  }

  {
    CGAL::Color_stream_guard warning(std::cout,
      std::vector<CGAL::Ansi_color>{CGAL::Ansi_color::Bold, CGAL::Ansi_color::Yellow});
    std::cout << "WARNING: Large dataset detected, may take longer\n";
  }

  {
    CGAL::Color_stream_guard processing(std::cout, CGAL::Ansi_color::Blue);
    std::cout << "Processing vertices...\n";
    std::cout << "Processing edges...\n";
    std::cout << "Processing faces...\n";
  }

  {
    CGAL::Color_stream_guard success(std::cout,
      std::vector<CGAL::Ansi_color>{CGAL::Ansi_color::Bold, CGAL::Ansi_color::Green});
    std::cout << "SUCCESS: Computation completed\n";
  }

  std::cout << "Algorithm finished\n\n";
}

void demonstrate_multiple_streams() {
  std::cout << "=== Multiple Streams with make_color_guards ===\n\n";

  std::cout << "Before coloring:\n";
  std::cout << "  cout: Normal output\n";
  std::cerr << "  cerr: Normal output\n";

  // Apply colors to both streams simultaneously
  {
    std::cout << "\nWith make_color_guards (Red):\n";
    auto guards = CGAL::make_color_guards(CGAL::Ansi_color::Red, std::cout, std::cerr);

    std::cout << "Red text on cout\n";
    std::cerr << "Red text on cerr\n";
  }

  std::cout << "\nAfter color scope:\n";
  std::cout << "  cout: Normal output again\n";
  std::cerr << "  cerr: Normal output again\n";

  // Using combined colors with multiple streams
  {
    std::cout << "\nBold green on multiple streams:\n";
    auto guards = CGAL::make_color_guards(
      std::vector<CGAL::Ansi_color>{CGAL::Ansi_color::Bold, CGAL::Ansi_color::Green},
      std::cout, std::clog);

    std::cout << "cout: Bold green output\n";
    std::clog << "clog: Bold green log message\n";
  }

  std::cout << "\n";
}

int main() {
  demonstrate_basic_colors();
  demonstrate_combined_colors();
  demonstrate_background_colors();
  demonstrate_nested_colors();
  simulate_colored_debug_output();
  demonstrate_multiple_streams();

  return 0;
}
