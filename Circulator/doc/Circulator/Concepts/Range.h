
/*!
\ingroup PkgHandlesAndCirculatorsConcepts
\cgalConcept

\cgal and the STL heavily use the concepts of iterators and iterator ranges
to describe linear sequences of elements, and algorithms operating on these.

The `Range` concept aims at encapsulating an iterator range, by providing
access to the first and past-the-end iterators of a range. The advantage is
that the syntax for passing ranges is much more concise than passing two
arguments separately.

Ranges come in different categories depending on the category of their iterator :
mutable or constant (modifiability of the elements pointed to), and forward,
bidirectional or random-access. The category can be queried using
`std::iterator_traits` and the corresponding iterator type. Note that
the concepts `Range` and `ConstRange` do not require anything on
the category or the value type of the iterator. It must be precised in the
documentation of any model of these concepts. For example, in the case of a vector of points, one would say:
<I>This type is a model of `Range` concept, its iterator type is random-access
and its value type is `Point`</I>.

Boost also offers the
<A HREF="https://www.boost.org/libs/range/">Boost.Range library</A>
which provides good support for ranges.

Finally, let us note that ranges, in general (especially in template context)
need to be passed and returned by (const) reference for efficiency. This is a
difference with iterators which are typically passed by value.

\cgalRefines `ConstRange`
\cgalRefines Boost's Range concept

\cgalHasModel STL containers

*/

class Range {
public:

/// \name Types
/// @{

/*!
The constant iterator type.
*/
typedef unspecified_type const_iterator;

/*!
The iterator type. It must be convertible to `const_iterator`.
*/
typedef unspecified_type iterator;

/*!
An unsigned integral type that can represent the
size of a range.
*/
typedef unspecified_type size_type;

/// @}

/// \name Member functions
/// @{

/*!
returns the const iterator pointing to the first element.
*/
const_iterator begin() const;

/*!
returns the past-the-end const iterator.
*/
const_iterator end() const;

/*!
returns the iterator pointing to the first element.
*/
iterator begin();

/*!
returns the past-the-end iterator.
*/
iterator end();

/*!
returns the size of the range.
*/
size_type size() const;

/*!
returns whether the range is empty.
*/
bool empty() const;

/// @}

}; /* end Range */

