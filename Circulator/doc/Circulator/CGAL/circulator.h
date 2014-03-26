namespace CGAL {


/*!
\addtogroup PkgHandlesAndCirculatorsAssert

Each of the following assertions, applicable to an iterator
or a circulator or both,
checks at compile time if its argument
is of the kind stated in the assertions name, i.e.\ a circulator, an
iterator, or a particular category of either an iterator or
a circulator. Note that neither input nor output circulators exists.

\sa `Circulator_tag`
\sa `Circulator_traits`
\sa `query_circulator_or_iterator`
\sa `Circulator`

 */


/*!
\ingroup PkgHandlesAndCirculatorsAssert

checks at compile time if its argument is a circulator.

*/
template <class C>
void Assert_circulator( const C &c);

/*!
\ingroup PkgHandlesAndCirculatorsAssert

checks at compile time if its argument is an iterator.

*/
template <class I>
void Assert_iterator( const I &i);

/*!
\ingroup PkgHandlesAndCirculatorsAssert

checks at compile time if its argument is a circulator or iterator.

*/
template< class IC>
void Assert_circulator_or_iterator(const IC& i);

/*!
\ingroup PkgHandlesAndCirculatorsAssert

*/
template <class I>
void Assert_input_category( const I &i);

/*!
\ingroup PkgHandlesAndCirculatorsAssert

*/
template <class I>
void Assert_output_category( const I &i);

/*!
\ingroup PkgHandlesAndCirculatorsAssert

*/
template <class IC>
void Assert_forward_category( const IC &ic);

/*!
\ingroup PkgHandlesAndCirculatorsAssert

*/
template <class IC>
void Assert_bidirectional_category( const IC &ic);

/*!
\ingroup PkgHandlesAndCirculatorsAssert

*/
template <class IC>
void Assert_random_access_category( const IC &ic);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgHandlesAndCirculatorsFunctions

The distance of a circulator `c` to a circulator `d` is the number of
elements in the range `[c, d)`. It is defined to be zero
for a circulator on an empty sequence and it returns the size of the data
structure when applied to a range of the form `[c, c)`.

\sa `circulator_size`
\sa `iterator_distance`
\sa `is_empty_range`
\sa `Circulator`

*/
template <class C> C::difference_type
circulator_distance(C c, C d);

} /* namespace CGAL */


namespace CGAL {

/*!
\ingroup PkgHandlesAndCirculatorsAdapter

The adaptor `Circulator_from_container` provides a circulator for an \stl container `C` of equal category as the iterator provided by the container.
The iterator must be at least of the forward iterator
category. The corresponding non-mutable circulator is called
`Const_circulator_from_container<C>`.

The container type `C` is supposed to conform to the \stl
requirements for container (i.e.\ to have a `begin()` and an
`end()` iterator as well as the local types
`reference`, `const_reference`, `value_type`,
`size_type`, and `difference_type`).

\cgalHeading{Types}

All types required for circulators are provided.

\cgalHeading{Operations}

The adaptor conforms to the requirements of the corresponding
circulator category. An additional member function
`current_iterator()` returns the current iterator pointing to
the same position as the circulator does.

\sa `Container_from_circulator`
\sa `Circulator_from_iterator`
\sa `Circulator`

\cgalHeading{Example}

The following program composes two adaptors - from a container to a
circulator and back to an iterator. It applies an \stl sort algorithm
on a \stl vector with three elements. The resulting vector will be
<TT>[2 5 9]</TT> as it is checked by the assertions. The program is
part of the \cgal distribution.

\cgalExample{Circulator/circulator_prog2.cpp}

*/
template< typename C >
class Circulator_from_container {
public:

/// \name Creation
/// @{

/*!
a circulator `c` on an empty sequence.
*/
Circulator_from_container();

/*!
a circulator `c` initialized to refer to the first element in
`container`, i.e.\ `container.begin()`.
The circulator `c` refers to an empty sequence if the
`container` is empty.

*/
Circulator_from_container(C* container);

/*!
a circulator `c` initialized to refer to the element `*i` in
`container`. \pre `*i` is dereferenceable and refers to `container`.

*/
Circulator_from_container(C* container, C::iterator i);

/// @}

}; /* end Circulator_from_container */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgHandlesAndCirculatorsAdapter

The adaptor `Circulator_from_iterator` converts two iterators of type
`I`, a begin and a past-the-end value, to a circulator of equal
category. The iterator must be at least of the forward iterator
category. The circulator will be mutable or non-mutable according to
the iterator. Iterators provide no `size_type`. This adapter
assumes `std::size_t` instead.

\cgalHeading{Operations}

The adaptor conforms to the requirements of the respective circulator
category. An additional member function `current_iterator()`
returns the current iterator pointing to the same position as the
circulator does.

\sa `Container_from_circulator`
\sa `Circulator_from_container`
\sa `Circulator`

\cgalHeading{Example}

The following program composes two adaptors - from an iterator to a
circulator and back to an iterator. It applies an \stl sort algorithm
on a \stl vector containing three elements. The resulting vector will
be <TT>[2 5 9]</TT> as it is checked by the assertions. The program is
part of the \cgal distribution.

\cgalExample{Circulator/circulator_prog1.cpp}

Another example usage for this adaptor is a random access circulator
over the built-in C arrays. Given an array of type <TT>T*</TT> with a
begin pointer <TT>b</TT> and a past-the-end pointer <TT>e</TT> the adaptor
`Circulator_from_iterator<T*> c(b,e)` is a random access circulator
`c` over this array.

*/
template< typename I >
class Circulator_from_iterator {
public:

/// \name Types
/// In addition all types required for circulators are provided.
/// @{

/*!

*/
typedef I iterator;

/// @}

/// \name Creation
/// @{

/*!
a circulator `c` on an empty sequence.
*/
Circulator_from_iterator();

/*!
a circulator `c` initialized to refer to the element
`*cur` in a range `[begin, end)`.
The circulator `c` refers to a empty sequence
if `begin==end`.

*/
Circulator_from_iterator(const I& begin,
const I& end, const I& cur = begin);

/*!
a copy of circulator `d` referring to the element `*cur`.
The circulator `c` refers to a empty sequence
if `d` does so.

*/
Circulator_from_iterator(
const Circulator_from_iterator<I,T,Size,Dist>& d,
const I& cur);

/// @}

}; /* end Circulator_from_iterator */
} /* end namespace CGAL */
namespace CGAL {

/*!
\ingroup PkgHandlesAndCirculatorsFunctions

The size of a circulator is the size of the data structure it refers
to. It is zero for a circulator on an empty sequence. The size can be
computed in linear time for forward and bidirectional circulators, and
in constant time for random access circulators using the minimal
circulator. The function `circulator_size(c)`
returns the circulator size. It uses the
`c.min_circulator()` function if `c` is a random
access circulator.

\sa `circulator_distance`
\sa `iterator_distance`
\sa `is_empty_range`
\sa `Circulator`

*/
template <class C> C::size_type circulator_size(C c);

} /* namespace CGAL */


namespace CGAL {

/*!
\ingroup PkgHandlesAndCirculators

The circulator traits class distinguishes between circulators and
iterators. It defines a local type `category` that is identical to the
type `Circulator_tag` if the iterator category of the argument
`C` is a circulator category. Otherwise it is identical to the type
`Iterator_tag`.

The local type `iterator_category` gives the corresponding
iterator category for circulators, i.e.\ one of
`forward_iterator_tag`, `bidirectional_iterator_tag`, or
`random_access_iterator_tag`.

The local type `circulator_category` gives the corresponding
circulator category for iterators, i.e.\ one of
`Forward_circulator_tag`, `Bidirectional_circulator_tag`, or
`Random_access_circulator_tag`.

\cgalHeading{Example}

A generic function `bar()` that distinguishes between a call with a
circulator range and a call with an iterator range:

\code{.cpp}

template <class I>
void bar( I i, I j, CGAL::Iterator_tag) {
  CGAL::Assert_iterator(i);
  // This function is called for iterator ranges [i,j).
}

template <class C>
void bar( C c, C d, CGAL::Circulator_tag) {
  CGAL::Assert_circulator(c);
  // This function is called for circulator ranges [c,d).
}

template <class IC>
void bar( IC i, IC j) { // calls the correct function
  return bar( i, j, typename CGAL::Circulator_traits<IC>::category());
}
\endcode

*/
template< typename C >
class Circulator_traits {
public:

/// \name Types
/// @{

/*!
either `Iterator_tag` or
`Circulator_tag`.
*/
typedef unspecified_type category;

/*!

corresponding iterator category for circulators.
*/
typedef unspecified_type iterator_category;

/*!

corresponding circulator category for iterator.
*/
typedef unspecified_type circulator_category;

/// @}

}; /* end Circulator_traits */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgHandlesAndCirculatorsAdapter

The adaptor `Container_from_circulator` is a class that converts any
circulator type `C` to a kind of container class, i.e.\ a class
that provides an `iterator` and a `const_iterator`
type and two member functions (`begin()` and `end()`) that return the appropriate iterators. By analogy to \stl container classes these member functions return a const iterator in
the case that the container itself is constant and a mutable iterator
otherwise.

\sa `Circulator_from_iterator`
\sa `Circulator_from_container`
\sa `Circulator`

\cgalHeading{Example}

The generic <TT>reverse()</TT> algorithm from the \stl can be used with an
adaptor if at least a bidirectional circulator <TT>c</TT> is given.

\code{.cpp}

Circulator c; // c is assumed to be a bidirectional circulator.
CGAL::Container_from_circulator<Circulator> container(c);
reverse( container.begin(), container.end());

\endcode

\cgalHeading{Implementation}

The iterator adaptor keeps track of the number of rounds a circulator
has done around the ring-like data structure (a kind of winding
number). It is used to distinguish between the start position and the
end position which will be denoted by the same circulator internally.
This winding number is zero for the `begin()`-iterator and one
for the `end()`-iterator. It is incremented whenever the
internal circulator passes the `begin()` position. Two
iterators are equal if their internally used circulators and winding
numbers are equal.
This is more general than necessary since an iterator equal to
`end()`-iterator is not supposed to be incremented
any more, which is here still possible in a defined manner.

The implementation is different for random access iterators.
The random access iterator has to be able to compute the size of the
data structure in constant time. This is for example needed if the
difference of the past-the-end iterator and the begin iterator is
taken, which is exactly the size of the data structure.
Therefore, if the circulator is of the random-access category, the
adapter chooses the minimal circulator for the internal anchor
position. The minimal circulator is part of the random access
circulator requirements, see
Page \ref sectionMinCircleRequ. For the random
access iterator the adaptor implements a total ordering relation that
is currently not required for random access circulators.

*/
template< typename C >
class Container_from_circulator {
public:

/// \name Types
/// @{

/*!

*/
typedef C Circulator;

/*!

*/
typedef unspecified_type iterator;

/*!

*/
typedef unspecified_type const_iterator;

/*!

*/
typedef unspecified_type value_type;

/*!

*/
typedef unspecified_type reference;

/*!

*/
typedef unspecified_type const_reference;

/*!

*/
typedef unspecified_type pointer;

/*!

*/
typedef unspecified_type const_pointer;

/*!

*/
typedef unspecified_type size_type;

/*!

*/
typedef unspecified_type difference_type;

/// @}

/// \name Creation
/// @{

/*!
any iterator of `container` will have a singular value.
*/
Container_from_circulator();

/*!
any iterator of `container` will have a singular value if the
circulator `c` corresponds to an empty sequence.
*/
Container_from_circulator(const C& c);

/// @}

/// \name Operations
/// The `iterator` and `const_iterator` types are of the appropriate
/// iterator category. In addition to the operations required for
/// their category, they have a member function `current_circulator()`
/// that returns a circulator pointing to the same position as the
/// iterator does.
/// @{

/*!
the start iterator.
*/
iterator begin();

/*!
the start const iterator.
*/
const_iterator begin() const;

/*!
the past-the-end iterator.
*/
iterator end();

/*!
the past-the-end const iterator.
*/
const_iterator end() const;

/// @}

}; /* end Container_from_circulator */
} /* end namespace CGAL */
namespace CGAL {

/*!
\ingroup PkgHandlesAndCirculatorsFunctions

is `true` if the range `[i, j)` is empty, `false` otherwise.

In order to write algorithms that work with iterator ranges as well as
with circulator ranges we have to consider the difference of
representing an empty range. For iterators this is the range `[i,i)`,
while for circulators it would be `c == NULL`, the empty sequence test.
The function `is_empty_range()` provides the necessary generic test
which accepts an iterator range or a circulator range and says whether
the range is empty or not.

\pre `IC` is either a circulator or an iterator type. The range `[i, j)` is valid.

\cgalHeading{Example}

The following function `process_all()` accepts a range `[i, j)` of an iterator or circulator `IC` and processes each
element in this range:

\code{.cpp}
template <class IC>
void process_all( IC i, IC j) {
  if (! CGAL::is_empty_range( i, j)) {
    do {
      process(*i);
    } while (++i != j);
  }
}
\endcode

\sa `iterator_distance`
\sa `CGAL_For_all`
\sa `Circulator_tag`
\sa `Circulator_traits`
\sa `Assert_circulator_or_iterator`
\sa `Circulator`

*/
template< class IC>
bool is_empty_range( const IC& i, const IC& j);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgHandlesAndCirculatorsFunctions

The following function returns the distance between either two
iterators or two circulators.

\sa `circulator_size`
\sa `circulator_distance`
\sa `is_empty_range`
\sa `Circulator_tag`
\sa `Assert_circulator_or_iterator`
\sa `CGAL_For_all`
\sa `Circulator`

*/
template <class IC> iterator_traits<IC>::difference_type
iterator_distance(IC ic1, IC ic2);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgHandlesAndCirculatorsFunctions

This function matches for type `I` if the iterator category of `I` belongs to an iterator.

\sa `Circulator_tag`
\sa `Circulator_traits`
\sa `Assert_circulator`
\sa `Circulator`
*/
template <class I>
Iterator_tag query_circulator_or_iterator( const I& i);

/*!
\ingroup PkgHandlesAndCirculatorsFunctions

This functiona matches for type `C` if the iterator category of `C` belongs to a circulator.

\sa `Circulator_tag`
\sa `Circulator_traits`
\sa `Assert_circulator`
\sa `Circulator`

*/
template <class C>
Circulator_tag query_circulator_or_iterator( const C& c);

} /* namespace CGAL */



/*!
\ingroup PkgHandlesAndCirculatorsFunctions

In order to write algorithms that work with iterator ranges as well as
with circulator ranges we have to consider the difference of
representing an empty range. For iterators this is the range `[i,i)`,
while for circulators it would be `c == NULL`, the empty sequence test.
The function `is_empty_range()` provides the necessary generic test
which accepts an iterator range or a circulator range and says whether
the range is empty or not.

A macro `CGAL_For_all( i, j)` simplifies the writing of such simple
loops as the one in the example of the function `is_empty_range()`.
`i` and `j` can be either iterators or circulators. The macro
loops through the range `[i, j)`. It increments `i` until it
reaches `j`. The implementation looks like:

\code
for ( bool _circ_loop_flag = ! ::CGAL::is_empty_range(i,j);
      _circ_loop_flag;
      _circ_loop_flag = ((++i) != (j))
)
\endcode


Note that the macro behaves like a `for`-loop. It can be used with
a single statement or with a statement block.  For bidirectional
iterators or circulators,  a backwards loop macro
`CGAL_For_all_backwards(i, j)` exists that decrements `j` until
it reaches `i`.

\sa `CGAL::iterator_distance()`
\sa `CGAL::is_empty_range()`
\sa `CGAL::Circulator_tag`
\sa `CGAL::Circulator_traits`
\sa `CGAL::Assert_circulator_or_iterator`
\sa `Circulator`

*/
#define CGAL_For_all(i,j)

/// \ingroup PkgHandlesAndCirculatorsFunctions
/// See ::CGAL_For_all
#define CGAL_For_all_backwards(i,j)
