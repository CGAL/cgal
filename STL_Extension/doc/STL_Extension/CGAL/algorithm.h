/// \defgroup STLAlgos Generic Algorithms
/// \ingroup PkgSTLExtensionRef


namespace CGAL {

/*!
\ingroup STLAlgos


Computes the minimal and the
maximal element of a range. It is modeled after the \stl functions
<a href="https://en.cppreference.com/w/cpp/algorithm/min_element">`std::min_element`</a>
and <a href="https://en.cppreference.com/w/cpp/algorithm/max_element">`std::max_element`</a>.
The advantage of `min_max_element()` compared to calling both \stl functions is that
one only iterates once over the sequence. This is more efficient
especially for large and/or complex sequences.

\cgalHeading{Example}

The following example program computes the minimal and
maximal element of the sequence `(3,\,6,\,5)`. Hence the output is
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
maximal element of a range. It is modeled after the \stl functions
<a href="https://en.cppreference.com/w/cpp/algorithm/min_element">`std::min_element`</a>
and <a
href="https://en.cppreference.com/w/cpp/algorithm/max_element">`std::max_element`</a>.
The advantage of `min_max_element()` compared to calling both \stl functions is that
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


namespace cpp98 {

/*!
\ingroup STLAlgos

Replacement for <a href="https://en.cppreference.com/w/cpp/algorithm/random_shuffle">`std::random_shuffle`</a>
which was deprecated in C++14, and removed by C++17.
In the \stl it was replaced by `std::shuffle`.

\note The implementation in \cgal produces the same order on all platforms.
*/
template <class RandomAccessIterator,
          class RandomGenerator>
void
random_shuffle(RandomAccessIterator begin, RandomAccessIterator end,
               RandomGenerator& random);
/*!
\ingroup STLAlgos

Replacement for <a href="https://en.cppreference.com/w/cpp/algorithm/random_shuffle">`std::random_shuffle`</a>
which was deprecated in C++14, and removed by C++17.
In the \stl it was replaced by `std::shuffle`.

\note The implementation in \cgal produces the same order on all platforms.
*/
template <class RandomAccessIterator>
void
random_shuffle(RandomAccessIterator begin, RandomAccessIterator end);
} // namespace cpp98

} /* namespace CGAL */
