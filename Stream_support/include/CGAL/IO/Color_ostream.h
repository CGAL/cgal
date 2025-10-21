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

namespace CGAL {

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
   * \param wrapped_buf The underlying streambuf to wrap
   * \param color The color to apply to the output
   */
  explicit
  Basic_color_streambuf(streambuf_type& wrapped_buf,
                        Ansi_color color = Ansi_color::Reset)
      : wrapped_buf_(&wrapped_buf)
      , color_code_(make_color_code(color))
      , at_line_start_(true)
  {
  }

  /**
   * \brief Construct a color streambuf wrapper with multiple colors.
   *
   * \param wrapped_buf The underlying streambuf to wrap
   * \param colors Vector of colors to combine (e.g., bold + red)
   */
  Basic_color_streambuf(streambuf_type& wrapped_buf,
                        const std::vector<Ansi_color>& colors)
      : wrapped_buf_(&wrapped_buf)
      , at_line_start_(true)
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
    }
    return *this;
  }

  ~Basic_color_streambuf() noexcept {
    if(wrapped_buf_) {
      if(!color_code_.empty() && !at_line_start_) {
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

    // If we're at the start of a line, output the color code first
    if(at_line_start_ && !color_code_.empty()) {
      if(!put_a_string(color_code_)) {
        return eof();
      }
      at_line_start_ = false;
    }

    // Check if this character is a newline
    if(to_char_type(ch) == char_type('\n')) {
      // Output reset code before newline
      if(!color_code_.empty()) {
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

} // namespace CGAL

#endif // CGAL_IO_COLOR_OSTREAM_H
