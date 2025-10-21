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

#ifndef CGAL_IO_INDENTING_OSTREAM_H
#define CGAL_IO_INDENTING_OSTREAM_H

#include <CGAL/config.h>

#include <ios>
#include <streambuf>
#include <string>
#include <tuple>

namespace CGAL {

/**
 * \ingroup PkgStreamSupportRef
 *
 * The class template `Basic_indenting_streambuf` wraps another `basic_streambuf`
 * and automatically adds indentation at the beginning of each line. This is
 * useful for formatting debug output with consistent indentation levels.
 *
 * \tparam CharT Character type (typically `char` or `wchar_t`)
 * \tparam Traits Character traits type
 */
template <typename CharT, typename Traits = std::char_traits<CharT>>
class Basic_indenting_streambuf : public std::basic_streambuf<CharT, Traits>, protected Traits
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
  string indent_string_;
  bool at_line_start_;

public:
  /**
   * \brief Construct an indenting streambuf wrapper.
   *
   * \param wrapped_buf The underlying streambuf to wrap
   * \param indent_string The string to use for indentation (default: 2 spaces)
   */
  explicit
  Basic_indenting_streambuf(streambuf_type& wrapped_buf,
                            const string& indent_string = {2, char_type(' ')})
      : wrapped_buf_(&wrapped_buf)
      , indent_string_(indent_string)
      , at_line_start_(true)
  {
  }

  /** \brief Get the current indentation string. */
  const string& indent_string() const { return indent_string_; }

  /**
   * \brief Set a new indentation string.
   *
   * \param new_indent The new indentation string
   */
  void set_indent_string(const string& new_indent) { indent_string_ = new_indent; }

  /**
   * \brief Set indentation level using repeated spaces.
   *
   * \param level Number of indentation levels
   * \param spaces_per_level Number of spaces per level (default: 2)
   */
  void set_indent_level(int level, int spaces_per_level = 2) {
    indent_string_ = string(level * spaces_per_level, char_type(' '));
  }

  /** \brief Get the wrapped streambuf. */
  streambuf_type& wrapped_streambuf() const { return *wrapped_buf_; }

protected:
  using traits_type::eof;
  using traits_type::not_eof;
  using traits_type::to_char_type;
  using traits_type::to_int_type;
  using traits_type::eq_int_type;

  int_type overflow(int_type ch = eof()) override {
    if(eq_int_type(ch, eof())) {
      return wrapped_buf_->pubsync() == 0 ? not_eof(ch) : eof();
    }

    // If we're at the start of a line, output the indentation first
    if(at_line_start_ && !indent_string_.empty()) {
      for(char_type indent_char : indent_string_) {
        if(eq_int_type(wrapped_buf_->sputc(indent_char), eof())) {
          return eof();
        }
      }
      at_line_start_ = false;
    }

    // Output the actual character
    int_type result = wrapped_buf_->sputc(to_char_type(ch));

    // Check if this character is a newline
    if(!eq_int_type(result, eof()) && to_char_type(ch) == char_type('\n')) {
      at_line_start_ = true;
    }

    return result;
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

/// Type alias for `Basic_indenting_streambuf<char>`
using Indenting_streambuf = Basic_indenting_streambuf<char>;

/// Type alias for `Basic_indenting_streambuf<wchar_t>`
using Indenting_wstreambuf = Basic_indenting_streambuf<wchar_t>;

/**
 * \ingroup PkgStreamSupportRef
 *
 * RAII helper class for temporarily installing an indenting streambuf on a stream.
 *
 * This class automatically restores the original streambuf when it goes out of scope.
 * It provides a convenient way to add indentation to output streams within a specific scope.
 *
 * \tparam StreamT The stream type (e.g., `std::ostream`, `std::wostream`)
 */
template <typename StreamT>
class Basic_indenting_stream_guard
{
private:
  using char_type = typename StreamT::char_type;
  using traits_type = typename StreamT::traits_type;
  using streambuf_type = Basic_indenting_streambuf<char_type, traits_type>;
  using wrapped_streambuf_type = std::basic_streambuf<char_type, traits_type>;
  using string = std::basic_string<char_type>;
  StreamT& stream_;
  wrapped_streambuf_type* original_buf_;
  string original_indent_string_;
  streambuf_type indenting_buf_;

public:
  /**
   * \brief Construct and install an indenting streambuf on the given stream.
   *
   * \param stream The stream to modify
   * \param indent_string The indentation string to use
   */
  explicit Basic_indenting_stream_guard(StreamT& stream,
                                        const string& indent_string)
      : stream_(stream)
      , original_buf_(stream.rdbuf())
      , indenting_buf_(*original_buf_, indent_string)
  {
    if(auto old_buf = dynamic_cast<streambuf_type*>(original_buf_)) {
      original_indent_string_ = old_buf->indent_string();
      old_buf->set_indent_string(original_indent_string_ + indent_string);
    } else {
      stream.rdbuf(&indenting_buf_);
    }
  }

  /**
   * \brief Construct and install an indenting streambuf on the given stream
   *
   * \param stream The stream to modify
   * \param spaces_per_level Number of indentation spaces
   */
  Basic_indenting_stream_guard(StreamT& stream, int spaces_per_level = 2)
      : Basic_indenting_stream_guard(stream,
                                     string(spaces_per_level, char_type(' '))) {
  }

  // Non-copyable, movable
  Basic_indenting_stream_guard(const Basic_indenting_stream_guard&) = delete;
  Basic_indenting_stream_guard& operator=(const Basic_indenting_stream_guard&) = delete;

  Basic_indenting_stream_guard(Basic_indenting_stream_guard&& other)
      : stream_(other.stream_)
      , original_buf_(other.original_buf_)
      , original_indent_string_(std::move(other.original_indent_string_))
      , indenting_buf_(std::move(other.indenting_buf_))
  {
    // Update the stream to point to the new buffer location if it was using the old one
    if(stream_.rdbuf() == &other.indenting_buf_) {
      stream_.rdbuf(&indenting_buf_);
    }
    // Mark other as moved-from so it won't restore the stream
    other.original_buf_ = nullptr;
  }

  Basic_indenting_stream_guard& operator=(Basic_indenting_stream_guard&&) = delete;

  /** \brief Destructor - restores the original streambuf. */
  ~Basic_indenting_stream_guard() {
    if(!original_buf_) return; // Moved-from object
    if(auto old_buf = dynamic_cast<streambuf_type*>(original_buf_)) {
      old_buf->set_indent_string(original_indent_string_);
    } else {
      stream_.rdbuf(original_buf_);
    }
  }
};

/// Type alias for `basic_indenting_stream_guard<std::ostream>`
using Indenting_stream_guard = Basic_indenting_stream_guard<std::ostream>;

/// Type alias for `basic_indenting_stream_guard<std::wostream>`
using Indenting_wstream_guard = Basic_indenting_stream_guard<std::wostream>;

/**
 * \ingroup PkgStreamSupportRef
 *
 * \brief Create indenting guards for multiple streams simultaneously.
 *
 * This helper function creates `Basic_indenting_stream_guard` objects for
 * multiple streams at once. All streams will use the same indentation settings,
 * and their original streambufs will be automatically restored when the guards
 * go out of scope. The guards are returned in a tuple.
 *
 * \tparam Streams The stream types (deduced from arguments)
 * \param spaces_per_level Number of spaces per indentation level
 * \param streams The streams to apply indentation to
 * \return A tuple of Basic_indenting_stream_guard objects
 *
 * Example usage:
 * \code
 * auto guards = CGAL::make_indenting_guards(2, std::cout, std::cerr);
 * std::cout << "Indented output\n";
 * std::cerr << "Indented error\n";
 * \endcode
 */
template <typename... Streams>
auto make_indenting_guards(int spaces_per_level, Streams&... streams) {
  return std::make_tuple(Basic_indenting_stream_guard<Streams>(streams, spaces_per_level)...);
}

/**
 * \ingroup PkgStreamSupportRef
 *
 * \brief Create indenting guards for multiple streams with a custom indent string.
 *
 * This overload allows specifying a custom indentation string instead of a number
 * of spaces.
 *
 * \tparam Streams The stream types (deduced from arguments)
 * \param indent_string The indentation string to use
 * \param streams The streams to apply indentation to
 * \return A tuple of Basic_indenting_stream_guard objects
 */
template <typename... Streams>
auto make_indenting_guards(const std::string& indent_string, Streams&... streams) {
  return std::make_tuple(Basic_indenting_stream_guard<Streams>(streams, indent_string)...);
}

} // namespace CGAL

#endif // CGAL_IO_INDENTING_OSTREAM_H