// Copyright (c) 2025  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_EXCEPTION_OSTREAM_H
#define CGAL_EXCEPTION_OSTREAM_H

#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace CGAL {

/**
 * \class Exception_basic_ostream
 * \brief A std::basic_ostream that throws an exception with its buffer content when flushed.
 *
 * Usage:
 * \code
 * CGAL::Exception_basic_ostream os;
 * os << "Error: " << value << std::endl; // throws std::runtime_error with the message
 * \endcode
 *
 * \tparam CharT Character type (default: char)
 * \tparam Traits Character traits (default: std::char_traits<CharT>)
 */
template <typename CharT = char, typename Traits = std::char_traits<CharT>>
class Exception_basic_ostream : public std::basic_ostream<CharT, Traits> {
  class buffer_type : public std::basic_stringbuf<CharT, Traits> {
  public:
    using int_type = typename Traits::int_type;
    buffer_type() = default;
    // When the buffer is flushed, throw an exception with the buffer content
    int sync() override {
      std::basic_string<CharT, Traits> msg = this->str();
      this->str({}); // clear buffer
      throw std::runtime_error(std::string(msg.begin(), msg.end()));
    }
  };
  buffer_type buffer_;
public:
  Exception_basic_ostream() : std::basic_ostream<CharT, Traits>(&buffer_) {}
  // Disallow copy and move
  Exception_basic_ostream(const Exception_basic_ostream&) = delete;
  Exception_basic_ostream& operator=(const Exception_basic_ostream&) = delete;
  Exception_basic_ostream(Exception_basic_ostream&&) = delete;
  Exception_basic_ostream& operator=(Exception_basic_ostream&&) = delete;
  ~Exception_basic_ostream() override = default;
};

/// /relates Exception_basic_ostream
using Exception_ostream = Exception_basic_ostream<char>;

/// /relates Exception_basic_ostream
using Exception_wostream = Exception_basic_ostream<wchar_t>;

inline Exception_ostream& exception_ostream() {
  static Exception_ostream os;
  return os;
}

} // namespace CGAL

#endif // CGAL_EXCEPTION_OSTREAM_H
