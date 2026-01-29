// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Efi Fogel       <efif@post.tau.ac.il>

#ifndef CGAL_VARIANT_OUTPUT_ITERATOR_H
#define CGAL_VARIANT_OUTPUT_ITERATOR_H

/* A compatible output iterator that accepts a std::variant and dispatches its
 * alternatives to different sinks via a visitor.
 *
 * Designed to replace boost::function_output_iterator when the algorithm emits
 * heterogeneous output (e.g., `Point_2` or `X_monotone_curve_2`).
 *
 * Many CGAL concepts (e.g., Make_x_monotone_2) require an operation that writes
 * results to an `OutputIterator` whose value type is:
 *
 * `std::variant<T1, T2, ...>`
 *
 * However:
 * 1. `boost::function_output_iterator` receives one type
 * 2. Visitors (operator()) are invoked after assignment
 * 3. CGAL freely copies, assigns, and default-constructs output iterators
 * `variant_output_iterator` bridges this gap.
 */

#include <variant>
#include <iterator>
#include <utility>
#include <memory>

namespace CGAL {
namespace detail {

template <class... Fs> struct overloaded : Fs... { using Fs::operator()...; };

template <class... Fs> overloaded(Fs...) -> overloaded<Fs...>;

} // namespace detail

template <typename Variant, typename Visitor>
class variant_output_iterator {
public:
  using value_type        = void;
  using difference_type   = std::ptrdiff_t;
  using iterator_category = std::output_iterator_tag;
  using pointer           = void;
  using reference         = void;

  /// Required by CGAL
  variant_output_iterator() = default;

  explicit variant_output_iterator(Visitor visitor) : visitor_(std::make_shared<Visitor>(std::move(visitor))) {}

  // Now ALWAYS copyable & assignable
  variant_output_iterator(const variant_output_iterator&) = default;
  variant_output_iterator& operator=(const variant_output_iterator&) = default;

  variant_output_iterator& operator*()     { return *this; }
  variant_output_iterator& operator++()    { return *this; }
  variant_output_iterator& operator++(int) { return *this; }

  variant_output_iterator& operator=(Variant const& v) {
    std::visit(*visitor_, v);
    return *this;
  }

private:
  std::shared_ptr<Visitor> visitor_;
};

template <typename Variant, typename... Fs>
auto make_variant_output_iterator(Fs&&... fs) {
  auto visitor = detail::overloaded{std::forward<Fs>(fs)...};
  return variant_output_iterator<Variant, decltype(visitor)>(std::move(visitor));
}

} // namespace CGAL

#endif
