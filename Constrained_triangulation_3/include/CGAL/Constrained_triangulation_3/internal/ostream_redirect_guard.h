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

#ifndef CGAL_CDT_3_INTERNAL_OSTREAM_REDIRECT_GUARD_H
#define CGAL_CDT_3_INTERNAL_OSTREAM_REDIRECT_GUARD_H

#include <ostream>

namespace CGAL {
namespace internal {

class ostream_redirect_guard {
public:
  // Builder for fluent API
  class builder {
  public:
    explicit builder(std::ostream& target) : target_(target) {}
    ostream_redirect_guard to(std::ostream& redirect_to) {
      return ostream_redirect_guard(target_, redirect_to);
    }
  private:
    std::ostream& target_;
  };

  static builder redirect(std::ostream& target) {
    return builder(target);
  }

  ostream_redirect_guard(const ostream_redirect_guard&) = delete;
  ostream_redirect_guard& operator=(const ostream_redirect_guard&) = delete;

  ~ostream_redirect_guard() { restore(); }

  ostream_redirect_guard(ostream_redirect_guard&& other) noexcept
    : target_(other.target_), old_buf_(other.old_buf_) {
    other.old_buf_ = nullptr;
  }
private:
  friend class builder;
  ostream_redirect_guard(std::ostream& target, std::ostream& redirect_to)
    : target_(target), old_buf_(target.rdbuf(redirect_to.rdbuf())) {}

  void restore() {
    if (old_buf_) target_.rdbuf(old_buf_);
    old_buf_ = nullptr;
  }

  std::ostream& target_;
  std::streambuf* old_buf_;
};

} // namespace internal
} // namespace CGAL

#endif // CGAL_CDT_3_INTERNAL_OSTREAM_REDIRECT_GUARD_H
