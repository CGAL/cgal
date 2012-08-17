
namespace CGAL {

/*!
\ingroup PkgHandlesAndCirculators

Iterators and circulators as well as different categories of
circulators can be distinguished with the use of discriminating
functions and the following circulator tags. A couple of base classes
simplify the task of writing own circulators. They declare the
appropriate tags and the local types needed for circulators.
To use the tags or base classes only it is sufficient to include:

\sa `query_circulator_or_iterator`
\sa `Circulator_traits`
\sa `Assert_circulator`
\sa `CGAL_For_all`
\sa `is_empty_range`
\sa `Circulator`

Example
--------------

The above declarations can be used to distinguish between iterators
and circulators and between different circulator categories. The
assertions can be used to protect a templatized algorithm against
instantiations that do not fulfill the requirements. The following
example program illustrates both.

\cgalexample{circulator_prog3.cpp}

Implementation
--------------

Since not all current compilers can eliminate the space needed for the
compile time tags even when deriving from them, we implement a variant
for each base class that contains a protected `void*` data member
called `_ptr`. Here, the allocated space in the derived
classes can be reused.

*/

class Circulator_tag {
public:

/// \name Compile Time Tags
/// @{

/*!
any circulator.
*/
struct Circulator_tag {};

/*!
any iterator.
*/
struct Iterator_tag {};

/*!

derived from `forward_iterator_tag`.
*/
struct Forward_circulator_tag {};

/*!

derived from `bidirectional_iterator_tag`.
*/
struct Bidirectional_circulator_tag {};

/*!

derived from `random_access_iterator_tag`.
*/
struct Random_access_circulator_tag {};

/// @}

/// \name Base Classes
/// @{

/*!

*/
template < class Category,
class T,
class Dist = std::ptrdiff_t,
class Size = std::size_t,
class Ptr = T*,
class Ref = T& >
struct Circulator_base {};

/*!

*/
template <class T, class Dist, class Size>
struct Forward_circulator_base {};

/*!

*/
template <class T, class Dist, class Size>
struct Bidirectional_circulator_base {};

/*!

*/
template <class T, class Dist, class Size>
struct Random_access_circulator_base {};

/// @}

/// \name Implementation
/// @{

/*!
forward circulator.
*/
template <class T, class Dist, class Size>
class Forward_circulator_ptrbase {};

/*!
bidirectional circulator.
*/
template <class T, class Dist, class Size>
class Bidirectional_circulator_ptrbase {};

/*!
random access circulator.
*/
template <class T, class Dist, class Size>
class Random_access_circulator_ptrbase {};

/// @}

}; /* end Circulator_tag */
} /* end namespace CGAL */
