/// \defgroup STLAlgos Generic Algorithms
/// \ingroup PkgStlExtension


namespace CGAL {

/*!
\ingroup STLAlgos

\deprecated This function is deprecated, CGAL::cpp11::copy_n should be
used instead.

Copies the first `n` items from `first` to `result`.

\returns the value of `result` after inserting the `n` items.

\note The \stl release June 13, 1997, from SGI contains an equivalent
function, but it is not part of the ISO standard.

\sa `CGAL::Counting_iterator<Iterator, Value>`

copies

*/
template <class InputIterator, class Size, class
OutputIterator> OutputIterator copy_n(InputIterator first, Size n,
OutputIterator result);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup STLAlgos


Computes the minimal and the
maximal element of a range. It is modeled after the STL functions
`std::min_element` and `std::max_element`. The advantage of
`min_max_element` compared to calling both STL functions is that
one only iterates once over the sequence. This is more efficient
especially for large and/or complex sequences.

\cgalHeading{Example}

The following example program computes the minimal and
maximal element of the sequence ` (3,\,6,\,5)`. Hence the output is
`min = 3, max = 6`.

\cgalExample{STL_Extension/min_max_element_example.cpp}

\returns a pair of iterators where
the first component refers to the minimal and the second component
refers to the maximal element in the range [`first`,
`last`). The ordering is defined by `operator<` on the
value type of `ForwardIterator`.
*/
template < class ForwardIterator > std::pair<
ForwardIterator, ForwardIterator > min_max_element(ForwardIterator
first, ForwardIterator last);


/*!
\ingroup STLAlgos

Computes the minimal and the
maximal element of a range. It is modeled after the STL functions
`std::min_element` and `std::max_element`. The advantage of
`min_max_element` compared to calling both STL functions is that
one only iterates once over the sequence. This is more efficient
especially for large and/or complex sequences.


\returns a pair of iterators where the first component refers to the minimal and the
second component refers to the maximal element in the range
[`first`, `last`).

\tparam CompareMin is an adaptable binary
function object: `VT` \f$ \times\f$ `VT` \f$ \rightarrow\f$ `bool` where `VT`
is the value type of `ForwardIterator`.

\tparam CompareMax is an adaptable binary
function object: `VT` \f$ \times\f$ `VT` \f$ \rightarrow\f$ `bool` where `VT`
is the value type of `ForwardIterator`.
*/
template < class ForwardIterator, class CompareMin,
class CompareMax > std::pair< ForwardIterator, ForwardIterator >
min_max_element(ForwardIterator first, ForwardIterator last,
CompareMin comp_min, CompareMax comp_max);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup STLAlgos

\deprecated This function is deprecated. `CGAL::cpp11::prev` should be used
instead.

Returns the previous iterator,
i.e.\ the result of `operator--` on a bidirectional iterator.

\sa `CGAL::successor()`

\returns `--it`.
*/
template <class BidirectionalIterator>
BidirectionalIterator predecessor(BidirectionalIterator it);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup STLAlgos

\deprecated This function is deprecated. `CGAL::cpp11::next` should be used
instead.


Returns the next iterator, i.e.
the result of `operator++` on a forward iterator.


\sa `CGAL::predecessor()`

\returns `++it`.
*/
template <class ForwardIterator>
ForwardIterator successor(ForwardIterator it);

namespace cpp11 {

/*!
\ingroup STLAlgos

The function returns the result of `operator++` on a
`ForwardIterator`. The exact behaviour is described in Paragraph 24.4.4
of the C++ standard draft
<a href="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2011/n3242.pdf">N3242</a>.

\note There is actually no function in namespace `CGAL::cpp11` with this
name, but a using declaration which imports a function from another
namespace. By order of priority: the one in namespace `std` is used
(provided by C++0x), if not found, then the one in namespace `boost`
is used.



\sa <a href="http://www.boost.org/doc/libs/1_46_1/libs/utility/utility.htm#functions_next_prior">boost::next</a>
\sa `CGAL::cpp11::prev()`

*/
template <typename ForwardIterator>
Iterator next(ForwardIterator it);

/*!
\ingroup STLAlgos

The function returns the result of `operator--` on
a `BidirectionalIterator`. The exact behaviour is described in
Paragraph 24.4.4 of the C++ standard draft
<a href="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2011/n3242.pdf">N3242</a>.

\note If C++0x is available the function `std::prev` is imported into
the namespace `CGAL::cpp11`, otherwise `CGAL::cpp11::prev` is declared with the
signature as given in Paragraph 24.4.4 of the ISO C++ Standard
and forwarded to `boost::prior`.
*/
template <typename BidirectionalIterator>
Iterator prev(BidirectionalIterator it);


/*!
\ingroup STLAlgos

Copies `n` items from an
input iterator to an output iterator. Its exact behaviour is defined
in Paragraph 25.3.1 of the C++ standard draft
<a href="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2011/n3242.pdf">N3242</a>.

\note This provides an implementation of the standard function
`copy_n` from the C++0x standard. If `copy_n` is available
in the `std::` namespace a using declaration is used, otherwise
an alternative implementation from \cgal is used.
*/

template< class InputIterator, class Size, class OutputIterator>
OutputIterator copy_n(InputIterator first, Size count, OutputIterator result);

} /* namespace cpp11 */
} /* namespace CGAL */

