namespace CGAL {
namespace cpp23 {

/*!
\ingroup PkgSTLExtensionUtilities

Replacement for `std::expected` that is added in C++23.

This class is a copy of [`tl::expected`](https://github.com/TartanLlama/expected),
with the namespace renamed to `CGAL::cpp23`. Its interface is that of
[`std::expected`](https://en.cppreference.com/w/cpp/utility/expected).
*/
template <typename T, typename E>
class expected {};

/*!
\ingroup PkgSTLExtensionUtilities

Replacement for `std::unexpected` that is added in C++23.
Wrapper for the error value of a `CGAL::cpp23::expected`.

\see CGAL::cpp23::expected<T, E>
*/
template <typename E>
class unexpected {};

/*!
\ingroup PkgSTLExtensionUtilities

Replacement for `std::bad_expected_access` that is added in C++23.

\see CGAL::cpp23::expected<T, E>
*/
template <typename E>
class bad_expected_access {};

} // namespace cpp23
} // namespace CGAL
