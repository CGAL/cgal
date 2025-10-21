#include <CGAL/IO/Color_ostream.h>
#include <iostream>

/*!
 * \example Stream_support/color_ostream.cpp
 *
 * This example demonstrates how to use CGAL::IO::Basic_color_streambuf
 * and CGAL::IO::Basic_color_stream_guard for formatting output with
 * ANSI colors.
 */

using CGAL::IO::Ansi_color;

void demonstrate_basic_colors() {
  std::cout << "=== Basic Color Usage ===\n\n";

  // Example 1: Single color
  {
    std::cout << "1. Single color (Red):\n";
    CGAL::IO::Color_stream_guard guard(std::cout, Ansi_color::Red);

    std::cout << "This text is red\n";
    std::cout << "All lines in this scope are red\n";
  }
  std::cout << "Back to normal color\n\n";

  // Example 2: Different colors
  {
    std::cout << "2. Different colors:\n";

    {
      CGAL::IO::Color_stream_guard green(std::cout, Ansi_color::Green);
      std::cout << "Green text\n";
    }

    {
      CGAL::IO::Color_stream_guard blue(std::cout, Ansi_color::Blue);
      std::cout << "Blue text\n";
    }

    {
      CGAL::IO::Color_stream_guard yellow(std::cout, Ansi_color::Yellow);
      std::cout << "Yellow text\n";
    }
  }
  std::cout << "Back to normal\n\n";

  // Example 3: Bright colors
  {
    std::cout << "3. Bright colors:\n";

    {
      CGAL::IO::Color_stream_guard bright_red(std::cout, Ansi_color::BrightRed);
      std::cout << "Bright red text\n";
    }

    {
      CGAL::IO::Color_stream_guard bright_green(std::cout, Ansi_color::BrightGreen);
      std::cout << "Bright green text\n";
    }

    {
      CGAL::IO::Color_stream_guard bright_cyan(std::cout, Ansi_color::BrightCyan);
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
    CGAL::IO::Color_stream_guard guard(std::cout,
      std::vector<Ansi_color>{Ansi_color::Bold, Ansi_color::Red});

    std::cout << "This text is bold and red\n";
  }
  std::cout << "Back to normal\n\n";

  // Example 2: Underlined text
  {
    std::cout << "2. Underlined colors:\n";
    CGAL::IO::Color_stream_guard guard(std::cout,
      std::vector<Ansi_color>{Ansi_color::Underline, Ansi_color::Blue});

    std::cout << "This text is underlined and blue\n";
  }
  std::cout << "Back to normal\n\n";

  // Example 3: Bold + underline + color
  {
    std::cout << "3. Multiple attributes:\n";
    CGAL::IO::Color_stream_guard guard(std::cout,
      std::vector<Ansi_color>{Ansi_color::Bold, Ansi_color::Underline, Ansi_color::Green});

    std::cout << "This text is bold, underlined, and green\n";
  }
  std::cout << "Back to normal\n\n";
}

void demonstrate_background_colors() {
  std::cout << "=== Background Colors ===\n\n";

  {
    std::cout << "1. Red background:\n";
    CGAL::IO::Color_stream_guard guard(std::cout, Ansi_color::BgRed);
    std::cout << "Text with red background\n";
  }
  std::cout << "Back to normal\n\n";

  {
    std::cout << "2. Foreground + background:\n";
    CGAL::IO::Color_stream_guard guard(std::cout,
      std::vector<Ansi_color>{Ansi_color::Yellow, Ansi_color::BgBlue});
    std::cout << "Yellow text on blue background\n";
  }
  std::cout << "Back to normal\n\n";

  {
    std::cout << "3. Bold white on bright green background:\n";
    CGAL::IO::Color_stream_guard guard(std::cout,
      std::vector<Ansi_color>{Ansi_color::Bold, Ansi_color::White, Ansi_color::BgBrightGreen});
    std::cout << "Bold white on bright green\n";
  }
  std::cout << "Back to normal\n\n";
}

void demonstrate_nested_colors() {
  std::cout << "=== Nested Color Scopes ===\n\n";

  {
    CGAL::IO::Color_stream_guard outer(std::cout, Ansi_color::Blue);
    std::cout << "Outer scope: blue text\n";

    {
      CGAL::IO::Color_stream_guard inner(std::cout, Ansi_color::Red);
      std::cout << "Inner scope: red text\n";

      {
        CGAL::IO::Color_stream_guard innermost(std::cout, Ansi_color::Green);
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
    CGAL::IO::Color_stream_guard info(std::cout, Ansi_color::Cyan);
    std::cout << "INFO: Initializing data structures\n";
  }

  {
    CGAL::IO::Color_stream_guard success(std::cout,
      std::vector<Ansi_color>{Ansi_color::Bold, Ansi_color::Green});
    std::cout << "SUCCESS: Data loaded successfully\n";
  }

  {
    CGAL::IO::Color_stream_guard warning(std::cout,
      std::vector<Ansi_color>{Ansi_color::Bold, Ansi_color::Yellow});
    std::cout << "WARNING: Large dataset detected, may take longer\n";
  }

  {
    CGAL::IO::Color_stream_guard processing(std::cout, Ansi_color::Blue);
    std::cout << "Processing vertices...\n";
    std::cout << "Processing edges...\n";
    std::cout << "Processing faces...\n";
  }

  {
    CGAL::IO::Color_stream_guard success(std::cout,
      std::vector<Ansi_color>{Ansi_color::Bold, Ansi_color::Green});
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
    auto guards = CGAL::IO::make_color_guards(Ansi_color::Red, std::cout, std::cerr);

    std::cout << "Red text on cout\n";
    std::cerr << "Red text on cerr\n";
  }

  std::cout << "\nAfter color scope:\n";
  std::cout << "  cout: Normal output again\n";
  std::cerr << "  cerr: Normal output again\n";

  // Using combined colors with multiple streams
  {
    std::cout << "\nBold green on multiple streams:\n";
    auto guards = CGAL::IO::make_color_guards(
      std::vector<Ansi_color>{Ansi_color::Bold, Ansi_color::Green},
      std::cout, std::clog);

    std::cout << "cout: Bold green output\n";
    std::clog << "clog: Bold green log message\n";
  }

  std::cout << "\n";
}

void test_color_detection() {
  std::cout << "=== Color Detection Tests ===\n\n";

  // Test automatic detection in streambuf
  std::cout << "1. Testing automatic color detection in Basic_color_streambuf:\n";
  {
    CGAL::IO::Color_streambuf color_buf(*std::cout.rdbuf(), Ansi_color::Green);
    std::cout << "   Colors automatically detected as: " << (color_buf.colors_enabled() ? "ENABLED" : "DISABLED")
              << "\n";

    // Temporarily install the buffer to show it works
    auto* old_buf = std::cout.rdbuf(&color_buf);
    std::cout << "   (This text should be green only if colors were auto-detected)\n";
    std::cout.rdbuf(old_buf);
  }
  std::cout << "\n";

  // Test stdout color support function
  std::cout << "2. Testing stdout_supports_color() function:\n";
  if(CGAL::IO::stdout_supports_color()) {
    std::cout << "   stdout supports colors: YES\n";
    {
      CGAL::IO::Color_stream_guard guard(std::cout, Ansi_color::Green);
      std::cout << "   (This text should be green if color is truly supported)\n";
    }
  } else {
    std::cout << "   stdout supports colors: NO\n";
    std::cout << "   (Colors are disabled - may be due to NO_COLOR env var,\n";
    std::cout << "    redirection to file, or unsupported terminal)\n";
  }
  std::cout << "\n";

  // Test stderr color support
  std::cout << "3. Testing stderr_supports_color():\n";
  if(CGAL::IO::stderr_supports_color()) {
    std::cout << "   stderr supports colors: YES\n";
    {
      CGAL::IO::Color_stream_guard guard(std::cerr, Ansi_color::Yellow);
      std::cerr << "   (This text should be yellow if color is truly supported)\n";
    }
  } else {
    std::cout << "   stderr supports colors: NO\n";
  }
  std::cout << "\n";

  // Test stream_supports_color() with different streams
  std::cout << "4. Testing stream_supports_color() for different streams:\n";
  std::cout << "   std::cout:" << (CGAL::IO::stream_supports_color(std::cout) ? "YES" : "NO") << "\n";
  std::cout << "   std::cerr:" << (CGAL::IO::stream_supports_color(std::cerr) ? "YES" : "NO") << "\n";
  std::cout << "   std::clog: " << (CGAL::IO::stream_supports_color(std::clog) ? "YES" : "NO") << "\n";
  std::cout << "\n";

  // Demonstrate automatic coloring (no manual check needed!)
  std::cout << "5. Automatic coloring example (no manual checking!):\n";
  {
    CGAL::IO::Color_stream_guard guard(std::cout,
                                   std::vector<Ansi_color>{Ansi_color::Bold, Ansi_color::Cyan});
    std::cout << "   This message is automatically colored only if terminal supports it!\n";
    std::cout << "   No need to check stream_supports_color() manually - it's automatic!\n";
  }
  std::cout << "\n";

  // Demonstrate smart error/warning coloring
  std::cout << "6. Smart colored output (automatically adapts to terminal capabilities):\n";

  // Info message - colors applied automatically!
  std::cout << "   ";
  {
    CGAL::IO::Color_stream_guard guard(std::cout, Ansi_color::Cyan);
    std::cout << "INFO:";
  }
  std::cout << " Normal information message\n";

  // Success message - colors applied automatically!
  std::cout << "   ";
  {
    CGAL::IO::Color_stream_guard guard(std::cout,
                                   std::vector<Ansi_color>{Ansi_color::Bold, Ansi_color::Green});
    std::cout << "SUCCESS:";
  }
  std::cout << " Operation completed successfully\n";

  // Warning message - colors applied automatically!
  std::cout << "   ";
  {
    CGAL::IO::Color_stream_guard guard(std::cout,
                                   std::vector<Ansi_color>{Ansi_color::Bold, Ansi_color::Yellow});
    std::cout << "WARNING:";
  }
  std::cout << " Potential issue detected\n";

  // Error message (using cerr) - colors applied automatically!
  std::cerr << "   ";
  {
    CGAL::IO::Color_stream_guard guard(std::cerr,
                                   std::vector<Ansi_color>{Ansi_color::Bold, Ansi_color::Red});
    std::cerr << "ERROR:";
  }
  std::cerr << " Critical error occurred\n";
  std::cout << "\n";

  // Environment hints
  std::cout << "7. Environment information:\n";
  const char* term = std::getenv("TERM");
  const char* no_color = std::getenv("NO_COLOR");

  std::cout << "   TERM variable: " << (term ? term : "(not set)") << "\n";
  std::cout << "   NO_COLOR variable: " << (no_color ? "SET (colors disabled)" : "(not set)") << "\n";

#ifdef _WIN32
  std::cout << "   Platform: Windows\n";
  const char* ansicon = std::getenv("ANSICON");
  std::cout << "   ANSICON variable: " << (ansicon ? ansicon : "(not set)") << "\n";
#else
  std::cout << "   Platform: POSIX (Linux/macOS/Unix)\n";
#endif

  std::cout << "\n";
  std::cout << "   Tip: To disable colors, set the NO_COLOR environment variable:\n";
  std::cout << "        export NO_COLOR=1\n";
  std::cout << "   Tip: When redirecting to a file, colors are automatically disabled:\n";
  std::cout << "        ./color_ostream > output.txt\n";
  std::cout << "\n";
}

int main() {
  // First, test color detection
  test_color_detection();

  // Then demonstrate color features
  demonstrate_basic_colors();
  demonstrate_combined_colors();
  demonstrate_background_colors();
  demonstrate_nested_colors();
  simulate_colored_debug_output();
  demonstrate_multiple_streams();

  return 0;
}
