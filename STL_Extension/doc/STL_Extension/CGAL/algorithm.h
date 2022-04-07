/// \defgroup STLAlgos Generic Algorithms
/// \ingroup PkgSTLExtensionRef


namespace CGAL {

/*!
\ingroup STLAlgos

\deprecated This function is deprecated, std::copy_n should be
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

\deprecated This function is deprecated. `std::prev` should be used
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

\deprecated This function is deprecated. `std::next` should be used
instead.


Returns the next iterator, i.e.
the result of `operator++` on a forward iterator.


\sa `CGAL::predecessor()`

\returns `++it`.
*/
template <class ForwardIterator>
ForwardIterator successor(ForwardIterator it);

namespace cpp98 {

/*!
\ingroup STLAlgos

Replacement for <a href="https://en.cppreference.com/w/cpp/algorithm/random_shuffle">`std::random_shuffle()`</a>
which was deprecated in C++14, and removed by C++17.
In the \stl it was replaced by `std::shuffle()`.

\note The implementation in \cgal produces the same order on all platforms.
*/
template <class RandomAccessIterator,
          class RandomGenerator>
void
random_shuffle(RandomAccessIterator begin, RandomAccessIterator end,
               RandomGenerator& random);
/*!
\ingroup STLAlgos

Replacement for <a href="https://en.cppreference.com/w/cpp/algorithm/random_shuffle">`std::random_shuffle()`</a>
which was deprecated in C++14, and removed by C++17.
In the \stl it was replaced by `std::shuffle()`.

\note The implementation in \cgal produces the same order on all platforms.
*/
template <class RandomAccessIterator>
void
random_shuffle(RandomAccessIterator begin, RandomAccessIterator end);
} // namespace cpp98

} /* namespace CGAL */

