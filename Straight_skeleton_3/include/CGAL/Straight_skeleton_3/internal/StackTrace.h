// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   util/StackTrace.h
 * author Gernot Walzl
 * date   2012-02-28
 */

#ifndef UTIL_STACKTRACE_H
#define UTIL_STACKTRACE_H

#ifdef __linux__
#include <execinfo.h>
#include <cxxabi.h>
#endif
#include <iostream>
#include <string>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {

class StackTrace {
public:
  virtual ~StackTrace();
  static std::string demangle(const std::string& symbol_mangled);
  static void print(std::ostream& os);
};

StackTrace::~StackTrace() {
    // intentionally does nothing
}

std::string StackTrace::demangle(const std::string& symbol_mangled)
{
  std::string result;
#ifdef __linux__
  size_t pos_begin = symbol_mangled.find('(');
  size_t pos_end = symbol_mangled.find('+');
  if (pos_begin != std::string::npos && pos_end != std::string::npos && pos_begin < pos_end) {
      std::string symbol = symbol_mangled.substr(pos_begin+1, pos_end-pos_begin-1);
      int status;
      char * realname = abi::__cxa_demangle(symbol.c_str(), 0, 0, &status);
      if (realname) {
        result.append(symbol_mangled.substr(0, pos_begin+1));
        result.append(realname);
        result.append(symbol_mangled.substr(pos_end));
        free(realname);
      } else {
        result.append(symbol_mangled);
      }
  } else {
      result.append(symbol_mangled);
  }
#endif
  return result;
}

void StackTrace::print(std::ostream& os) {
#ifdef __linux__
  const int size = 100;
  void* buffer[size];
  int nptrs = backtrace(buffer, size);
  char** strings = backtrace_symbols(buffer, nptrs);
  for (int i = 1; i < nptrs-2; ++i) {
    os << demangle(strings[i]) << std::endl;
  }
  free(strings);
#endif
}

} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* UTIL_STACKTRACE_H */
