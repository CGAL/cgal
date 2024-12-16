namespace CGAL {

/*!
\ingroup PkgSTLExtensionRef

\returns `std::array<T, N>` where `N` is the number of arguments given to
the function. The position of each argument in the array is the same
as its position in the argument list.
*/
template <class T> std::array<T, N> make_array(const T&...);

/*!

Functor that constructs `std::array<T, N>` where `N` is the number of
arguments given to the function. The position of each argument in the
array is the same as its position in the argument list.

This is the functor version of `make_array()`.
*/
struct Construct_array
{
  template <class T> std::array<T, N> operator()(const T&...);
};

} /* end namespace CGAL */
