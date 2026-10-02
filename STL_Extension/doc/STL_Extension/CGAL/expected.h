// Copyright (c) 2025 GeometryFactory (France)
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : The CGAL Project

namespace CGAL {
namespace cpp23 {

/*!
\ingroup PkgSTLExtensionRef

\brief A type representing either an expected value or an unexpected error.

An object of type `CGAL::cpp23::expected<T, E>` holds either:
- a *value* of type `T` (the expected case), or
- an *error* of type `E` wrapped in `CGAL::cpp23::unexpected<E>` (the unexpected case).

This is \cgal's vendored implementation of `std::expected` as standardized
in C++23 (ISO/IEC 14882:2024). It is provided as a compatibility shim for
compilers that do not yet fully support C++23's `<expected>` header.
The implementation is sourced from the
<a href="https://github.com/TartanLlama/expected">TartanLlama/expected</a> library
(v1.3.1, CC0 license) with the following changes:
- namespace renamed from `tl` to `CGAL::cpp23`,
- internal macros prefixed with `CGAL_TL_` to avoid collisions.

Include `<CGAL/expected.h>` to use this type.

\sa `CGAL::cpp23::unexpected<E>`
\sa `CGAL::cpp23::bad_expected_access<E>`
\sa <a href="https://en.cppreference.com/w/cpp/utility/expected"><code>std::expected</code> on cppreference</a>

\cgalHeading{Example}

The following shows how to use `CGAL::cpp23::expected` to return a
result-or-error from a parsing function, avoiding exceptions:

\code{.cpp}
#include <CGAL/expected.h>
#include <string>
#include <charconv>

// Returns the integer value, or an error string on failure.
CGAL::cpp23::expected<int, std::string>
parse_int(std::string_view s)
{
  int value{};
  auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
  if (ec != std::errc{})
    return CGAL::cpp23::unexpected<std::string>("parse error: " + std::string(s));
  return value;
}

int main()
{
  auto r = parse_int("42");
  if (r.has_value())
    std::cout << "parsed: " << r.value() << "\n";
  else
    std::cout << "error: " << r.error() << "\n";
}
\endcode

\tparam T The type of the expected value. May be `void` to represent a
          computation that either succeeds (producing no value) or fails with
          an error.
\tparam E The error type. Must not be `void`.
*/
template <class T, class E>
class expected
{
public:
  /// \name Types
  /// @{
  typedef T value_type;  ///< The type of the expected value.
  typedef E error_type;  ///< The type of the error.
  typedef unexpected<E> unexpected_type;  ///< Convenience alias.
  /// @}

  /// \name Constructors
  /// @{
  /*!
  \brief Default constructor. Constructs an object containing a
         default-constructed value of type `T`.
  \pre `T` is default-constructible.
  */
  constexpr expected();

  /*!
  \brief Constructs an object containing an error from an `unexpected<E>` object.
  \param u An unexpected wrapper holding the error value.
  */
  template <class G>
  constexpr expected(const unexpected<G>& u);

  /*!
  \brief Constructs an object containing an unexpected error value in-place,
         using the `unexpect` tag.
  \param tag `CGAL::cpp23::unexpect`
  \param args Arguments forwarded to the constructor of `E`.
  */
  template <class... Args>
  constexpr explicit expected(unexpect_t tag, Args&&... args);
  /// @}

  /// \name Observers
  /// @{
  /*!
  Returns `true` if `*this` holds an expected value, `false` if it holds an error.
  */
  constexpr bool has_value() const noexcept;

  /*!
  Equivalent to `has_value()`.
  */
  constexpr explicit operator bool() const noexcept;

  /*!
  Returns a reference to the contained expected value.
  If `has_value()` is `false`, this function throws
  `CGAL::cpp23::bad_expected_access<E>` (when exceptions are enabled),
  or calls `std::terminate()` otherwise.
  \throws CGAL::cpp23::bad_expected_access<E> if `has_value()` is `false`
          and exceptions are enabled.
  */
  constexpr const T& value() const &;

  /*!
  Returns a reference to the contained error.
  \pre `has_value()` is `false`.
  */
  constexpr const E& error() const &;

  /*!
  Returns the contained value if `has_value()`, otherwise returns `default_value`.
  \tparam U A type convertible to `T`.
  */
  template <class U>
  constexpr T value_or(U&& default_value) const &;
  /// @}

  /// \name Monadic Operations
  /// @{
  /*!
  If `has_value()`, calls `f(value())` and returns its result (an `expected`);
  otherwise propagates the error.
  */
  template <class F>
  constexpr auto and_then(F&& f) &;

  /*!
  If `has_value()`, applies `f` to the contained value and wraps the result in
  an `expected`; otherwise propagates the error unchanged.
  */
  template <class F>
  constexpr auto map(F&& f) &;

  /*!
  If `!has_value()`, calls `f(error())` and returns its result;
  otherwise propagates the value.
  */
  template <class F>
  constexpr auto or_else(F&& f) &;

  /*!
  If `!has_value()`, applies `f` to the error and returns an `expected` with
  the transformed error; otherwise propagates the value unchanged.
  */
  template <class F>
  constexpr auto map_error(F&& f) &;

  /*!
  C++23 standard-name alias for `map()`. If `has_value()`, applies `f` to the
  contained value and wraps the result in an `expected`; otherwise propagates
  the error unchanged.
  */
  template <class F>
  constexpr auto transform(F&& f) &;

  /*!
  C++23 standard-name alias for `map_error()`. If `!has_value()`, applies `f`
  to the error and returns an `expected` with the transformed error; otherwise
  propagates the value unchanged.
  */
  template <class F>
  constexpr auto transform_error(F&& f) &;
  /// @}
};

/*!
\ingroup PkgSTLExtensionRef

\brief Wraps an error value for use with `CGAL::cpp23::expected`.

An `unexpected<E>` holds an error of type `E`. It is used to
construct an `expected<T,E>` object that is in the error state:
\code{.cpp}
auto err = CGAL::cpp23::unexpected<std::string>("something went wrong");
CGAL::cpp23::expected<int, std::string> r = err;
assert(!r.has_value());
\endcode

\tparam E The error type. Must not be `void`.

\sa `CGAL::cpp23::expected<T,E>`
\sa <a href="https://en.cppreference.com/w/cpp/utility/expected/unexpected"><code>std::unexpected</code> on cppreference</a>
*/
template <class E>
class unexpected
{
public:
  /*!
  Constructs an `unexpected` wrapping a copy of `e`.
  */
  constexpr explicit unexpected(const E& e);

  /*!
  Constructs an `unexpected` by moving from `e`.
  */
  constexpr explicit unexpected(E&& e);

  /*!
  Returns a reference to the stored error value.
  */
  constexpr const E& value() const &;
};

/*!
\ingroup PkgSTLExtensionRef

\brief Exception thrown when calling `expected<T,E>::value()` on an object
       in the error state.

Thrown by `CGAL::cpp23::expected<T,E>::value()` when the object does not hold
an expected value (i.e., `has_value()` is `false`). Only thrown when
exceptions are enabled (`__EXCEPTIONS` or `_CPPUNWIND` is defined);
otherwise `std::terminate()` is called.

\sa `CGAL::cpp23::expected<T,E>`
\sa <a href="https://en.cppreference.com/w/cpp/utility/expected/bad_expected_access"><code>std::bad_expected_access</code> on cppreference</a>
*/
template <class E>
class bad_expected_access : public std::exception
{
public:
  /*!
  Constructs a `bad_expected_access` wrapping error value `e`.
  \param e The error value that caused the bad access.
  */
  explicit bad_expected_access(E e);

  /*!
  Returns a human-readable description of the error.
  Always returns `"Bad expected access"`.
  */
  virtual const char* what() const noexcept override;

  /*!
  Returns a reference to the stored error value.
  */
  constexpr const E& error() const &;
};

/*!
\ingroup PkgSTLExtensionRef

\brief Tag type used to construct `expected<T,E>` in the error state in-place.

\sa `CGAL::cpp23::unexpect`
*/
struct unexpect_t {};

/*!
\ingroup PkgSTLExtensionRef

\brief Tag constant of type `unexpect_t`.

Pass this constant as the first argument to the `expected<T,E>` constructor
to construct an error value in-place:
\code{.cpp}
CGAL::cpp23::expected<int, std::string> r(CGAL::cpp23::unexpect, "error msg");
assert(!r.has_value());
\endcode

\sa `CGAL::cpp23::unexpect_t`
*/
constexpr unexpect_t unexpect;

/*!
\ingroup PkgSTLExtensionRef

\brief Creates an `unexpected<E>` from a value.

Convenience factory equivalent to `unexpected<std::decay_t<E>>(std::forward<E>(e))`.

\tparam E The error type (deduced).
\param   e The error value.
\returns An `unexpected<std::decay_t<E>>` wrapping `e`.

\sa `CGAL::cpp23::unexpected<E>`
*/
template <class E>
unexpected<typename std::decay<E>::type> make_unexpected(E&& e);

} /* end namespace cpp23 */
} /* end namespace CGAL */
