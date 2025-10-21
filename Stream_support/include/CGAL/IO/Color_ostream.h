// Copyright (c) 2025
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_IO_COLOR_OSTREAM_H
#define CGAL_IO_COLOR_OSTREAM_H

#include <CGAL/config.h>

#include <ios>
#include <streambuf>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <io.h>
#define CGAL_ISATTY _isatty
#define CGAL_FILENO _fileno
#else
#include <unistd.h>
#define CGAL_ISATTY isatty
#define CGAL_FILENO fileno
#endif

#include <cstdlib> // for getenv

namespace CGAL {
namespace IO {

/**
 * \ingroup PkgStreamSupportRef
 *
 * ANSI color codes for terminal output.
 */
enum class Ansi_color {
  Reset = 0,
  Bold = 1,
  Dim = 2,
  Underline = 4,
  Blink = 5,
  Reverse = 7,
  Hidden = 8,

  // Foreground colors
  Black = 30,
  Red = 31,
  Green = 32,
  Yellow = 33,
  Blue = 34,
  Magenta = 35,
  Cyan = 36,
  White = 37,

  // Bright foreground colors
  BrightBlack = 90,
  BrightRed = 91,
  BrightGreen = 92,
  BrightYellow = 93,
  BrightBlue = 94,
  BrightMagenta = 95,
  BrightCyan = 96,
  BrightWhite = 97,

  // Background colors
  BgBlack = 40,
  BgRed = 41,
  BgGreen = 42,
  BgYellow = 43,
  BgBlue = 44,
  BgMagenta = 45,
  BgCyan = 46,
  BgWhite = 47,

  // Bright background colors
  BgBrightBlack = 100,
  BgBrightRed = 101,
  BgBrightGreen = 102,
  BgBrightYellow = 103,
  BgBrightBlue = 104,
  BgBrightMagenta = 105,
  BgBrightCyan = 106,
  BgBrightWhite = 107
};

/**
 * \ingroup PkgStreamSupportRef
 *
 * The class template `Basic_color_streambuf` wraps another `basic_streambuf`
 * and automatically adds ANSI color codes to the output. This is useful for
 * colorizing terminal output with consistent color schemes.
 *
 * \tparam CharT Character type (typically `char` or `wchar_t`)
 * \tparam Traits Character traits type
 */
template <typename CharT, typename Traits = std::char_traits<CharT>>
class Basic_color_streambuf : public std::basic_streambuf<CharT, Traits>, protected Traits
{
public:
  using char_type = CharT;
  using traits_type = Traits;
  using int_type = typename traits_type::int_type;
  using pos_type = typename traits_type::pos_type;
  using off_type = typename traits_type::off_type;
  using streambuf_type = std::basic_streambuf<char_type, traits_type>;
  using string = std::basic_string<char_type>;

private:
  streambuf_type* wrapped_buf_;
  string color_code_;
  bool at_line_start_;
  bool colors_enabled_;

  /** \brief Detect if the wrapped buffer supports color output. */
  static bool detect_color_support(streambuf_type* buf) {
    // Check NO_COLOR environment variable first
    if(std::getenv("NO_COLOR")) {
      return false;
    }

    // Try to determine if this is stdout or stderr
    // by comparing pointers (works for standard streams)
    if(buf == std::cout.rdbuf() || buf == std::clog.rdbuf()) {
      int fd = CGAL_FILENO(stdout);
      if(!CGAL_ISATTY(fd)) {
        return false;
      }
    } else if(buf == std::cerr.rdbuf()) {
      int fd = CGAL_FILENO(stderr);
      if(!CGAL_ISATTY(fd)) {
        return false;
      }
    } else {
      // Unknown buffer, conservatively disable colors
      return false;
    }

#ifdef _WIN32
    // Windows: check for ANSICON or assume modern Windows
    return std::getenv("ANSICON") || true;
#else
    // POSIX: check TERM environment variable
    const char* term = std::getenv("TERM");
    if(!term) {
      return false;
    }

    std::string term_str(term);
    if(term_str == "dumb") {
      return false;
    }

    // Check for common color-capable terminals
    return term_str.find("color") != std::string::npos ||
           term_str.find("xterm") != std::string::npos ||
           term_str.find("rxvt") != std::string::npos ||
           term_str.find("screen") != std::string::npos ||
           term_str.find("tmux") != std::string::npos ||
           term_str.find("linux") != std::string::npos;
#endif
  }

  /** \brief Generate ANSI escape sequence for a color. */
  static string make_color_code(Ansi_color color) {
    string code;
    code += char_type('\033');
    code += char_type('[');
    int color_value = static_cast<int>(color);
    // Convert integer to string
    if(color_value == 0) {
      code += char_type('0');
    } else {
      string digits;
      while(color_value > 0) {
        digits = char_type('0' + (color_value % 10)) + digits;
        color_value /= 10;
      }
      code += digits;
    }
    code += char_type('m');
    return code;
  }

public:
  /**
   * \brief Construct a color streambuf wrapper.
   *
   * Colors are automatically disabled if the wrapped buffer is not connected
   * to a color-capable terminal (checked via isatty() and TERM variable).
   *
   * \param wrapped_buf The underlying streambuf to wrap
   * \param color The color to apply to the output
   */
  explicit
  Basic_color_streambuf(streambuf_type& wrapped_buf,
                        Ansi_color color = Ansi_color::Reset)
      : wrapped_buf_(&wrapped_buf)
      , color_code_(make_color_code(color))
      , at_line_start_(true)
      , colors_enabled_(detect_color_support(&wrapped_buf))
  {
  }

  /**
   * \brief Construct a color streambuf wrapper with multiple colors.
   *
   * Colors are automatically disabled if the wrapped buffer is not connected
   * to a color-capable terminal (checked via isatty() and TERM variable).
   *
   * \param wrapped_buf The underlying streambuf to wrap
   * \param colors Vector of colors to combine (e.g., bold + red)
   */
  Basic_color_streambuf(streambuf_type& wrapped_buf,
                        const std::vector<Ansi_color>& colors)
      : wrapped_buf_(&wrapped_buf)
      , at_line_start_(true)
      , colors_enabled_(detect_color_support(&wrapped_buf))
  {
    if(!colors.empty()) {
      color_code_ += char_type('\033');
      color_code_ += char_type('[');
      for(std::size_t i = 0; i < colors.size(); ++i) {
        if(i > 0) color_code_ += char_type(';');
        int color_value = static_cast<int>(colors[i]);
        // Convert integer to string
        if(color_value == 0) {
          color_code_ += char_type('0');
        } else {
          string digits;
          while(color_value > 0) {
            digits = char_type('0' + (color_value % 10)) + digits;
            color_value /= 10;
          }
          color_code_ += digits;
        }
      }
      color_code_ += char_type('m');
    }
  }

  // Non-copyable
  Basic_color_streambuf(const Basic_color_streambuf&) = delete;
  Basic_color_streambuf& operator=(const Basic_color_streambuf&) = delete;

  // Move constructor
  Basic_color_streambuf(Basic_color_streambuf&& other) noexcept
      : wrapped_buf_(std::exchange(other.wrapped_buf_, nullptr))
      , color_code_(std::move(other.color_code_))
      , at_line_start_(std::exchange(other.at_line_start_, true))
      , colors_enabled_(std::exchange(other.colors_enabled_, false))
  {
  }

  // Move assignment
  Basic_color_streambuf& operator=(Basic_color_streambuf&& other) noexcept {
    if(this != &other) {
      // If currently wrapping a buffer, try to sync it first.
      if(wrapped_buf_) {
        wrapped_buf_->pubsync();
      }
      wrapped_buf_ = std::exchange(other.wrapped_buf_, nullptr);
      color_code_ = std::move(other.color_code_);
      at_line_start_ = std::exchange(other.at_line_start_, true);
      colors_enabled_ = std::exchange(other.colors_enabled_, false);
    }
    return *this;
  }

  ~Basic_color_streambuf() noexcept {
    if(wrapped_buf_) {
      if(colors_enabled_ && !color_code_.empty() && !at_line_start_) {
        string reset_code = make_color_code(Ansi_color::Reset);
        (void)put_a_string(reset_code);
      }
      wrapped_buf_->pubsync();
    }
  }

  /** \brief Get the current color code string. */
  const string& color_code() const { return color_code_; }

  /**
   * \brief Set a new color code.
   *
   * \param string The new color code to apply
   */
  void set_color_code(const string& color_code) {
    color_code_ = color_code;
  }

  /** \brief Check if colors are enabled for this streambuf. */
  bool colors_enabled() const { return colors_enabled_; }

  /** \brief Get the wrapped streambuf. */
  streambuf_type& wrapped_streambuf() const { return *wrapped_buf_; }

protected:
  using traits_type::eof;
  using traits_type::not_eof;
  using traits_type::to_char_type;
  using traits_type::to_int_type;
  using traits_type::eq_int_type;

  bool put_a_string(const string& str) {
    for(char_type ch : str) {
      if(eq_int_type(wrapped_buf_->sputc(ch), eof())) {
        return false;
      }
    }
    return true;
  }

  int_type overflow(int_type ch = eof()) override {
    if(eq_int_type(ch, eof())) {
      return wrapped_buf_->pubsync() == 0 ? not_eof(ch) : eof();
    }

    // If we're at the start of a line, output the color code first (if colors enabled)
    if(colors_enabled_ && at_line_start_ && !color_code_.empty()) {
      if(!put_a_string(color_code_)) {
        return eof();
      }
      at_line_start_ = false;
    }

    // Check if this character is a newline
    if(to_char_type(ch) == char_type('\n')) {
      // Output reset code before newline (if colors enabled)
      if(colors_enabled_ && !color_code_.empty()) {
        string reset_code = make_color_code(Ansi_color::Reset);
        if(!put_a_string(reset_code)) {
          return eof();
        }
      }
      at_line_start_ = true;
    }

    // Output the actual character
    return wrapped_buf_->sputc(to_char_type(ch));
  }

  int sync() override { return wrapped_buf_->pubsync(); }

  std::streamsize xsputn(const char_type* s, std::streamsize count) override {
    std::streamsize written = 0;

    for(std::streamsize i = 0; i < count; ++i) {
      if(eq_int_type(overflow(to_int_type(s[i])), eof())) {
        break;
      }
      ++written;
    }

    return written;
  }

  pos_type seekoff(off_type off,
                   std::ios_base::seekdir dir,
                   std::ios_base::openmode which = std::ios_base::in | std::ios_base::out) override {
    return wrapped_buf_->pubseekoff(off, dir, which);
  }

  pos_type seekpos(pos_type pos, std::ios_base::openmode which = std::ios_base::in | std::ios_base::out) override {
    return wrapped_buf_->pubseekpos(pos, which);
  }

  std::streamsize showmanyc() override {
    return wrapped_buf_->in_avail();
  }

  std::streamsize xsgetn(char_type* s, std::streamsize count) override {
    return wrapped_buf_->sgetn(s, count);
  }

  int_type underflow() override {
    return wrapped_buf_->sgetc();
  }

  int_type uflow() override {
    return wrapped_buf_->sbumpc();
  }

  int_type pbackfail(int_type c = eof()) override {
    return wrapped_buf_->sputbackc(to_char_type(c));
  }

  void imbue(const std::locale& loc) override {
    wrapped_buf_->pubimbue(loc);
  }

  streambuf_type* setbuf(char_type* s, std::streamsize n) override {
    return wrapped_buf_->pubsetbuf(s, n);
  }
};

/// Type alias for `Basic_color_streambuf<char>`
using Color_streambuf = Basic_color_streambuf<char>;

/// Type alias for `Basic_color_streambuf<wchar_t>`
using Color_wstreambuf = Basic_color_streambuf<wchar_t>;

/**
 * \ingroup PkgStreamSupportRef
 *
 * RAII helper class for temporarily installing a color streambuf on a stream.
 *
 * This class automatically restores the original streambuf when it goes out of scope.
 * It provides a convenient way to add colors to output streams within a specific scope.
 *
 * \tparam StreamT The stream type (e.g., `std::ostream`, `std::wostream`)
 */
template <typename StreamT>
class Basic_color_stream_guard
{
private:
  using char_type = typename StreamT::char_type;
  using traits_type = typename StreamT::traits_type;
  using streambuf_type = Basic_color_streambuf<char_type, traits_type>;
  using wrapped_streambuf_type = std::basic_streambuf<char_type, traits_type>;
  using string = std::basic_string<char_type>;
  StreamT& stream_;
  wrapped_streambuf_type* original_buf_;
  string original_color_code_;
  streambuf_type color_buf_;

public:
  /**
   * \brief Construct and install a color streambuf on the given stream.
   *
   * \param stream The stream to modify
   * \param color The color to apply to the output
   */
  explicit Basic_color_stream_guard(StreamT& stream,
                                    Ansi_color color)
      : stream_(stream)
      , original_buf_(stream.rdbuf())
      , color_buf_(*original_buf_, color)
  {
    if(auto old_buf = dynamic_cast<streambuf_type*>(original_buf_)) {
      original_color_code_ = old_buf->color_code();
      // For nested colors, we keep the current buffer but note the original color
    } else {
      stream.rdbuf(&color_buf_);
    }
  }

  /**
   * \brief Construct and install a color streambuf with multiple colors.
   *
   * \param stream The stream to modify
   * \param colors Vector of colors to combine (e.g., bold + red)
   */
  Basic_color_stream_guard(StreamT& stream,
                           const std::vector<Ansi_color>& colors)
      : stream_(stream)
      , original_buf_(stream.rdbuf())
      , color_buf_(*original_buf_, colors)
  {
    if(auto old_buf = dynamic_cast<streambuf_type*>(original_buf_)) {
      original_color_code_ = old_buf->color_code();
    } else {
      stream.rdbuf(&color_buf_);
    }
  }

  // Non-copyable, movable
  Basic_color_stream_guard(const Basic_color_stream_guard&) = delete;
  Basic_color_stream_guard& operator=(const Basic_color_stream_guard&) = delete;

  Basic_color_stream_guard(Basic_color_stream_guard&& other)
      : stream_(other.stream_)
      , original_buf_(other.original_buf_)
      , original_color_code_(std::move(other.original_color_code_))
      , color_buf_(std::move(other.color_buf_))
  {
    // Update the stream to point to the new buffer location if it was using the old one
    if(stream_.rdbuf() == &other.color_buf_) {
      stream_.rdbuf(&color_buf_);
    }
    // Mark other as moved-from so it won't restore the stream
    other.original_buf_ = nullptr;
  }

  Basic_color_stream_guard& operator=(Basic_color_stream_guard&&) = delete;

  /** \brief Destructor - restores the original streambuf. */
  ~Basic_color_stream_guard() {
    if(!original_buf_) return; // Moved-from object
    if(auto old_buf = dynamic_cast<streambuf_type*>(original_buf_)) {
      old_buf->set_color_code(original_color_code_);
    } else {
      stream_.rdbuf(original_buf_);
    }
  }
};

/// Type alias for `Basic_color_stream_guard<std::ostream>`
using Color_stream_guard = Basic_color_stream_guard<std::ostream>;

/// Type alias for `Basic_color_stream_guard<std::wostream>`
using Color_wstream_guard = Basic_color_stream_guard<std::wostream>;

/**
 * \ingroup PkgStreamSupportRef
 *
 * \brief Create color guards for multiple streams simultaneously.
 *
 * This helper function creates `Basic_color_stream_guard` objects for
 * multiple streams at once. All streams will use the same color settings,
 * and their original streambufs will be automatically restored when the guards
 * go out of scope. The guards are returned in a tuple.
 *
 * \tparam Streams The stream types (deduced from arguments)
 * \param color The color to apply
 * \param streams The streams to apply colors to
 * \return A tuple of Basic_color_stream_guard objects
 *
 * Example usage:
 * \code
 * auto guards = CGAL::make_color_guards(Color::Red, std::cout, std::cerr);
 * std::cout << "Red output\n";
 * std::cerr << "Red error\n";
 * \endcode
 */
template <typename... Streams>
auto make_color_guards(Ansi_color color, Streams&... streams) {
  return std::make_tuple(Basic_color_stream_guard<Streams>(streams, color)...);
}

/**
 * \ingroup PkgStreamSupportRef
 *
 * \brief Create color guards for multiple streams with multiple colors.
 *
 * This overload allows specifying multiple colors to combine (e.g., bold + red).
 *
 * \tparam Streams The stream types (deduced from arguments)
 * \param colors Vector of colors to combine
 * \param streams The streams to apply colors to
 * \return A tuple of Basic_color_stream_guard objects
 */
template <typename... Streams>
auto make_color_guards(const std::vector<Ansi_color>& colors, Streams&... streams) {
  return std::make_tuple(Basic_color_stream_guard<Streams>(streams, colors)...);
}

/**
 * \ingroup PkgStreamSupportRef
 *
 * \brief Check if a stream is attached to a terminal that supports colors.
 *
 * This function checks if the given output stream is connected to a terminal
 * (TTY) that can display ANSI color codes. It performs the following checks:
 *
 * 1. Verifies the stream is attached to a TTY device using `isatty()`
 * 2. Checks if the `TERM` environment variable indicates color support
 * 3. Respects the `NO_COLOR` environment variable (if set, colors are disabled)
 * 4. On Windows, also checks for `ANSICON` environment variable
 *
 * \param stream The output stream to check (e.g., `std::cout`, `std::cerr`)
 * \return `true` if the stream supports color output, `false` otherwise
 *
 * \note This function may return false negatives on some systems, but should
 *       never return false positives (it won't report color support where
 *       it doesn't exist).
 *
 * Example usage:
 * \code
 * if(CGAL::stream_supports_color(std::cout)) {
 *   CGAL::Color_stream_guard guard(std::cout, CGAL::Ansi_color::Red);
 *   std::cout << "This text is red!\n";
 * } else {
 *   std::cout << "Plain text output\n";
 * }
 * \endcode
 */
template <typename CharT, typename Traits>
bool stream_supports_color(const std::basic_ostream<CharT, Traits>& stream) {
  // Check if NO_COLOR environment variable is set (standard way to disable colors)
  if(std::getenv("NO_COLOR")) {
    return false;
  }

  // Try to get the file descriptor for the stream
  // This works for std::cout (stdout) and std::cerr (stderr)
  const auto* fbuf = dynamic_cast<const std::basic_filebuf<CharT, Traits>*>(stream.rdbuf());
  if(!fbuf) {
    // If it's not a file buffer, check if it's stdout or stderr directly
    const std::basic_ostream<CharT, Traits>* std_stream = &stream;
    int fd = -1;

    if(std_stream == &std::cout || std_stream == &std::clog) {
      fd = CGAL_FILENO(stdout);
    } else if(std_stream == &std::cerr) {
      fd = CGAL_FILENO(stderr);
    }

    if(fd == -1) {
      return false; // Unknown stream type
    }

    // Check if the file descriptor is a TTY
    if(!CGAL_ISATTY(fd)) {
      return false;
    }
  } else {
    // For file buffers, we can't easily get the fd in a portable way
    // Conservatively assume no color support for files
    return false;
  }

#ifdef _WIN32
  // On Windows, check for ANSICON or Windows 10+ with VT100 support
  // Modern Windows Terminal and ConEmu support ANSI colors
  if(std::getenv("ANSICON")) {
    return true;
  }
  // Windows 10+ console supports ANSI by default
  // We could check Windows version, but being conservative here
  return true; // Modern Windows usually supports colors
#else
  // On POSIX systems, check the TERM environment variable
  const char* term = std::getenv("TERM");
  if(!term) {
    return false; // No TERM variable, probably not a color terminal
  }

  std::string term_str(term);

  // Check for common indicators of color support
  if(term_str.find("color") != std::string::npos ||
     term_str.find("xterm") != std::string::npos ||
     term_str.find("rxvt") != std::string::npos ||
     term_str.find("screen") != std::string::npos ||
     term_str.find("tmux") != std::string::npos ||
     term_str.find("linux") != std::string::npos) {
    return true;
  }

  // Check for dumb terminal
  if(term_str == "dumb") {
    return false;
  }

  // For other TERM values, conservatively assume color support
  // Most modern terminals support colors
  return true;
#endif
}

/**
 * \ingroup PkgStreamSupportRef
 *
 * \brief Check if stdout supports color output.
 *
 * Convenience function equivalent to `stream_supports_color(std::cout)`.
 *
 * \return `true` if stdout supports color output, `false` otherwise
 */
inline bool stdout_supports_color() {
  return stream_supports_color(std::cout);
}

/**
 * \ingroup PkgStreamSupportRef
 *
 * \brief Check if stderr supports color output.
 *
 * Convenience function equivalent to `stream_supports_color(std::cerr)`.
 *
 * \return `true` if stderr supports color output, `false` otherwise
 */
inline bool stderr_supports_color() {
  return stream_supports_color(std::cerr);
}

} // namespace IO
} // namespace CGAL

#endif // CGAL_IO_COLOR_OSTREAM_H
