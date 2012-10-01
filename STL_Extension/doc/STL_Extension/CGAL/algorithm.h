/// \defgroup STLAlgos Generic Algorithms
/// \ingroup PkgStlExtension


namespace CGAL {

/*!
\ingroup STLAlgos

\deprecated This function is deprecated, CGAL::cpp0x::copy_n should be
used instead.

The function `copy_n` copies \f$ n\f$ items from an input iterator to
an output iterator which is useful for possibly infinite sequences of
random geometric objects.

\note The \stl release June 13, 1997, from SGI contains an equivalent
function, but it is not part of the ISO standard.



\sa `CGAL::Counting_iterator<Iterator, Value>` 

copies the first \f$ n\f$ items from `first` to `result`. Returns the
value of `result` after inserting the \f$ n\f$ items.

*/
template <class InputIterator, class Size, class
OutputIterator> OutputIterator copy_n(InputIterator first, Size n,
OutputIterator result);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup STLAlgos


The function `min_max_element` computes the minimal and the 
maximal element of a range. It is modeled after the STL functions 
`min_element` and `max_element`. The advantage of 
`min_max_element` compared to calling both STL functions is that 
one only iterates once over the sequence. This is more efficient 
especially for large and/or complex sequences. 

### Example ###

The following example program computes the minimal and 
maximal element of the sequence \f$ (3,\,6,\,5)\f$. Hence the output is 
`min = 3, max = 6`. 

\cgalexample{STL_Extension/min_max_element_example.cpp} 

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

The function `min_max_element` computes the minimal and the 
maximal element of a range. It is modeled after the STL functions 
`min_element` and `max_element`. The advantage of 
`min_max_element` compared to calling both STL functions is that 
one only iterates once over the sequence. This is more efficient 
especially for large and/or complex sequences. 

### Example ###

The following example program computes the minimal and 
maximal element of the sequence \f$ (3,\,6,\,5)\f$. Hence the output is 
`min = 3, max = 6`. 

\cgalexample{STL_Extension/min_max_element_example.cpp} 

\returns a pair of iterators where the first component refers to the minimal and the
second component refers to the maximal element in the range
[`first`, `last`). 

\requires `CompareMin` and `CompareMax` are adaptable binary
function objects: `VT` \f$ \times\f$ `VT` \f$ \rightarrow\f$ `bool` where `VT`
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

\deprecated This function is deprecated. `CGAL::cpp0x::prev` should be used 
instead. 

The function `predecessor` returns the previous iterator, 
i.e. the result of `operator--` on a bidirectional iterator. 

\sa `CGAL::successor` 

\returns `--it`.
*/
template <class BidirectionalIterator>
BidirectionalIterator predecessor(BidirectionalIterator it);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup STLAlgos

\deprecated This function is deprecated. `CGAL::cpp0x::next` should be used 
instead. 

The function `successor` returns the next iterator, i.e. 
the result of `operator++` on a forward iterator. 



\sa `CGAL::predecessor` 

\returns `++it`.
*/
template <class ForwardIterator>
ForwardIterator successor(ForwardIterator it);

namespace cpp0x {

/*!
\ingroup STLAlgos

The function returns the result of `operator++` on a
ForwardIterator. The exact behaviour is described in \f$ \mathsection 24.4.4 \f$ 
of the C++ standard draft 
<a href="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2011/n3242.pdf">N3242</a>.

\note There is actually no function in namespace `CGAL::cpp0x` with this
name, but a using declaration which imports a function from another
namespace. By order of priority: the one in namespace `std` is used
(provided by C++0x), if not found, then the one in namespace `boost`
is used.



\sa <a href="http://www.boost.org/doc/libs/1_46_1/libs/utility/utility.htm#functions_next_prior">boost::next</a>
\sa CGAL::cpp0x::prev

*/
template <typename ForwardIterator>
Iterator next(ForwardIterator it);

/*!
\ingroup STLAlgos

The function returns the result of `operator--` on
a BidirectionalIterator. The exact behaviour is described in
\f$\mathsection 24.4.4\f$ of the C++ standard draft
<a href="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2011/n3242.pdf">N3242</a>.

\note If C++0x is available the function `std::prev` is imported into
the namespace `CGAL::cpp0x`, otherwise `CGAL::cpp0x::prev` is declared with the
signature as given in \f$\mathsection 24.4.4\f$ of the ISO C++ Standard
and forwarded to `boost::prior`.
*/
template <typename BidirectionalIterator>
Iterator prev(BidirectionalIterator it);


/*!
\ingroup STLAlgos

The function `copy_n` copies `n` items from an
input iterator to an output iterator. Its exact behaviour is defined
in \f$\mathsection 25.3.1\f$ of the C++ standard draft
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

