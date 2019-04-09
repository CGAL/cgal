
namespace CGAL {

/*!
\ingroup PkgSTLExtensionUtilities


The Quadruple class is an extension of
`std::pair`. `Quadruple` is a heterogeneous quadruple: it holds
one object of type `T1`, one of type `T2`, one of type
`T3`, and one of type `T4`. A `Quadruple` is much like a
container, in that it "owns" its elements. It is not actually a
model of container, though, because it does not support the standard
methods (such as iterators) for accessing the elements of a
container.

`std::tuple` or `std::array` instead for new uses.



\cgalHeading{Requirements}

`T1`, `T2`, `T3` and `T4` must be
`Assignable`. Additional operations have additional requirements.

*/
template< typename T1, typename T2, typename T3, typename T4 >
class Quadruple {
public:

/// \name Types
/// @{
/*!

*/
typedef T1 first_type;



/// @}


/// \name Types
/// @{
/*!

*/
typedef T2 second_type;



/// @}


/// \name Types
/// @{
/*!

*/
typedef T3 third_type;



/// @}


/// \name Types
/// @{
/*!

*/
typedef T4 fourth_type;



/// @}


/// \name Variables
/// @{
/*!
first element. Please access it using `get<0>()`.
*/
T1 first;



/// @}


/// \name Variables
/// @{
/*!
second element. Please access it using `get<1>()`.
*/
T2 second;



/// @}


/// \name Variables
/// @{
/*!
third element. Please access it using `get<2>()`.
*/
T3 third;



/// @}


/// \name Variables
/// @{
/*!
fourth element. Please access it using `get<3>()`.
*/
T4 fourth;



/// @}


/// \name Creation
/// @{
/*!
introduces a quadruple using the
default constructor of the four elements.
*/
Quadruple();



/// @}


/// \name Creation
/// @{
/*!
constructs a
quadruple such that `first` is constructed from `x`,
`second` is constructed from `y`, `third` is
constructed from `z`, and `fourth` is constructed from
`w`.
*/
Quadruple(T1 x, T2 y, T3 z, T4 w);



/// @}


/// \name Creation
/// @{
/*!
constructs a quadruple such that
`first` is constructed from `u`, `second` is
constructed from `v`, `third` is constructed from `w`,
and `fourth` is constructed from `x`.

Proper conversion operators must exist from `U` to `T1`, `V` to `T2`, `W` to `T3`, and `X` to `T4`.
*/
template <class U, class V, class W, class X>
Quadruple(U u, V v, W w, X x);



/// @}


/// \name Creation
/// @{
/*!
Gives access to `first`, `second`, `third` or `fourth`
whenever `i` is 0, 1, 2 or 3, via a, potentially const, reference.
Note: `T` stands for `T1`, `T2`, `T3` or `T4`
depending on `i`.
*/
template <int i> T get();



/// @}

/*!
The comparison operator. It uses lexicographic comparison:
the return value is true if the first element of `x` is less
than the first element of `y`, and false if the first element
of `y` is less than the first element of `x`. If neither
of these is the case, then it returns true if the second element
of `x` is less than the second element of `y`, and false
if the second element of `y` is less than the second element
of `x`. If neither of these is the case, then it returns true
if the third element of `x` is less than the third element of
`y`, and false if the third element of `y` is less than
the third element of `x`. If neither of these is the case,
then it returns the result of comparing the fourth elements of
`x` and `y`. This operator may only be used if `T1`,
`T2`, `T3`, and `T4` define the comparison operator.
*/
template <class T1, class T2, class T3, class T4> bool
operator<(Quadruple<T1, T2, T3, T4> x, Quadruple<T1, T2, T3, T4>
y);

/*!
The equality operator. The return value is true if and only
the first elements of `x` and `y` are equal, the second
elements of `x` and `y` are equal, the third elements of
`x` and `y` are equal, and the fourth elements of `x`
and `y` are equal. This operator may only be used if
`T1`, `T2`, `T3`, and `T4` define the equality
operator.
*/
template <class T1, class T2, class T3, class T4> bool
operator==(Quadruple<T1, T2, T3, T4> x, Quadruple<T1, T2, T3, T4>
y);


}; /* end Quadruple */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSTLExtensionUtilities



The Triple class is an extension of `std::pair`.
`Triple` is a heterogeneous triple: it holds one object of type
`T1`, one of type `T2`, and one of type `T3`. A
`Triple` is much like a container, in that it "owns" its
elements. It is not actually a model of container, though, because
it does not support the standard methods (such as iterators) for
accessing the elements of a container.

`std::tuple` or `std::array` instead for new uses.



\cgalHeading{Requirements}

`T1`, `T2` and `T3` must be `Assignable`.
Additional operations have additional requirements.

*/
template< typename T1, typename T2, typename T3 >
class Triple {
public:

/// \name Types
/// @{
/*!

*/
typedef T1 first_type;



/// @}


/// \name Types
/// @{
/*!

*/
typedef T2 second_type;



/// @}


/// \name Types
/// @{
/*!

*/
typedef T3 third_type;



/// @}


/// \name Variables
/// @{
/*!
first element. Please access it using `get<0>()`.
*/
T1 first;



/// @}


/// \name Variables
/// @{
/*!
second element. Please access it using `get<1>()`.
*/
T2 second;



/// @}


/// \name Variables
/// @{
/*!
third element. Please access it using `get<2>()`.
*/
T3 third;



/// @}


/// \name Creation
/// @{
/*!
introduces a triple using the default
constructor of the three elements.
*/
Triple();



/// @}


/// \name Creation
/// @{
/*!
constructs a triple such
that `first` is constructed from `x`, `second` is
constructed from `y`, and `third` is constructed from
`z`.
*/
Triple(T1 x, T2 y, T3 z);



/// @}


/// \name Creation
/// @{
/*!
constructs a triple such that `first` is constructed
from `u`, `second` is constructed from `v`, and
`third` is constructed from `w`.
Proper conversion operators must exist from `U` to `T1`, `V` to `T2`, and `W` to `T3`.
*/
template <class U, class V, class W> Triple(U u, V v,
W w);



/// @}


/// \name Creation
/// @{
/*!
Gives access to `first`, `second` or `third` whenever
`i` is 0, 1 or 2, via a, potentially const, reference.
Note: `T` stands for `T1`, `T2` or `T3` depending
on `i`.
*/
template <int i> T get();



/// @}

/*!
The
comparison operator. It uses lexicographic comparison: the return
value is true if the first element of `x` is less than the
first element of `y`, and false if the first element of
`y` is less than the first element of `x`. If neither of
these is the case, then it returns true if the second element of
`x` is less than the second element of `y`, and false if
the second element of `y` is less than the second element of
`x`. If neither of these is the case, then it returns the
result of comparing the third elements of `x` and `y`.
This operator may only be used if `T1`, `T2` and `T3`
define the comparison operator.
*/
template <class T1, class T2, class T3> bool
operator<(Triple<T1, T2, T3> x, Triple<T1, T2, T3> y);

/*!

The
equality operator. The return value is true if and only the first
elements of `x` and `y` are equal, the second elements of
`x` and `y` are equal, and the third elements of `x`
and `y` are equal. This operator may only be used if
`T1`, `T2` and `T3` define the equality operator.
*/
template <class T1, class T2, class T3> bool
operator==(Triple<T1, T2, T3> x, Triple<T1, T2, T3> y);

}; /* end Triple */

/*!
Equivalent to `Quadruple<T1, T2, T3, T4>(x, y, z, w)`.
\relates Quadruple
*/
template <class T1, class T2, class T3, class T4>
Quadruple<T1, T2, T3, T4> make_quadruple(T1 x, T2 y, T3 z, T4 w);

/*!
Equivalent to `Quadruple<T1, T2, T3, T4>(x, y, z, w)`.
\relates Quadruple
*/
template <class T1, class T2, class T3, class T4>
Quadruple<T1, T2, T3, T4> make_tuple(T1 x, T2 y, T3 z, T4 w);

#ifndef DOXYGEN_RUNNING
/*!
\ingroup PkgSTLExtensionUtilities
Creates a pair `(t1,t2)` if `comp(t1,t2)==true` and `(t2,t1)` otherwise.
`comp` is a binary function taking two elements of type T that returns
a value convertible to `bool`.
An overload is provided with `std::less<T>` as default for `Compare`.
\note If rvalue references are supported by the compiler, the return type of
the function is that of `std::make_pair(t1,t2)`.
*/
template <class T, class Compare>
std::pair<T,T> make_sorted_pair(T t1, T t2, Compare comp)
{
  return comp(t1, t2) ? std::make_pair(t1,t2) : std::make_pair(t2,t1);
}
#endif

} /* end namespace CGAL */
