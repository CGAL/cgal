
/*!
\ingroup PkgHandlesAndCirculatorsConcepts
\cgalConcept

Most data structures in \cgal use the concept of `Handle` in their user
interface to refer to the elements they store. This concept describes what is
sometimes called a trivial iterator. A `Handle` is akin to a pointer to
an object providing the dereference operator `operator*()` and member
access `operator->()` but no increment or decrement operators like
iterators. A `Handle` is intended to be used whenever the referenced
object is not part of a logical sequence.

Like iterators, the handle can be passed as template argument to
`std::iterators_traits` in order to extract its `value_type`,
the type of the element pointed to.
The `iterator_category` is `void`.

\cgalRefines{Descriptor}

The default constructed object must be unique as far as the equality
operator is concerned (this serves the same purpose as NULL for pointers).
(Note that this is not a generally supported feature of iterators of
standard containers.)

\cgalHasModelsBegin
\cgalHasModels{T* (pointer)}
\cgalHasModels{const T* (const pointers)}
\cgalHasModels{Iterator}
\cgalHasModels{Circulator}
\cgalHasModelsEnd
*/
class Handle {
public:

/// \name Dereference
/// @{

/*!
returns the object pointed to.
*/
value_type& operator*();

/*!
returns a pointer to the object pointed to.
*/
value_type* operator->();

/// @}

}; /* end Handle */

