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

#include <ios>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>

#include <CGAL/exceptions.h>

namespace CGAL {

/**
 * \class Exception_basic_ostream
 * \brief A stream-like object that throws an exception with its buffer content when destroyed.
 *
 * Usage:
 * \code
 * {
 *   CGAL::Exception_basic_ostream os;
 *   os << "Error: " << value;
 * } // throws std::runtime_error with the message when os goes out of scope
 * \endcode
 *
 * \note This class is move-only.
 *
 * \tparam CharT Character type (default: `char`)
 * \tparam Traits Character traits (default: `std::char_traits<CharT>`)
 */
template <typename CharT = char, typename Traits = std::char_traits<CharT>>
class Exception_basic_ostream {
  std::basic_ostringstream<CharT, Traits> stream_;
public:
  Exception_basic_ostream() : stream_() {
    stream_.precision(std::cerr.precision());
  }

  // move-only
  Exception_basic_ostream(Exception_basic_ostream&&) = default;
  Exception_basic_ostream& operator=(Exception_basic_ostream&&) = default;

  ~Exception_basic_ostream() noexcept(false) {
    std::basic_string<CharT, Traits> msg = stream_.str();
    if(!msg.empty()) {
      throw CGAL::Failure_exception("CGAL", "", __FILE__, __LINE__, std::string(msg.begin(), msg.end()));
    }
  }

  template<typename T>
  Exception_basic_ostream& operator<<(T&& value) {
    stream_ << std::forward<T>(value);
    return *this;
  }

  // Support for stream manipulators
  Exception_basic_ostream& operator<<(std::basic_ostream<CharT, Traits>& (*manip)(std::basic_ostream<CharT, Traits>&)) {
    stream_ << manip;
    return *this;
  }

  Exception_basic_ostream& operator<<(std::basic_ios<CharT, Traits>& (*manip)(std::basic_ios<CharT, Traits>&)) {
    stream_ << manip;
    return *this;
  }

  Exception_basic_ostream& operator<<(std::ios_base& (*manip)(std::ios_base&)) {
    stream_ << manip;
    return *this;
  }

};

/// /relates Exception_basic_ostream
using Exception_ostream = Exception_basic_ostream<char>;

/// /relates Exception_basic_ostream
using Exception_wostream = Exception_basic_ostream<wchar_t>;

inline Exception_ostream exception_ostream() {
  Exception_ostream os;
  return os;
}

} // namespace CGAL

#endif // CGAL_EXCEPTION_OSTREAM_H
