
namespace CGAL {

/*!
\addtogroup PkgHandlesAndCirculatorsTags

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

\cgalHeading{Example}

The above declarations can be used to distinguish between iterators
and circulators and between different circulator categories. The
assertions can be used to protect a templatized algorithm against
instantiations that do not fulfill the requirements. The following
example program illustrates both.

\cgalExample{Circulator/circulator_prog3.cpp}

*/


/*!
\addtogroup PkgHandlesAndCirculatorsBaseClasses

\cgalHeading{Implementation}

Since not all current compilers can eliminate the space needed for the
compile time tags even when deriving from them, we implement a variant
for each base class that contains a protected `void*` data member
called `_ptr`. Here, the allocated space in the derived
classes can be reused.


*/


/*!
\ingroup PkgHandlesAndCirculatorsTags
A tag for any circulator type.
*/
struct Circulator_tag {};

/*!
\ingroup PkgHandlesAndCirculatorsTags
A tag for any iterator type.
*/
struct Iterator_tag {};

/*!
\ingroup PkgHandlesAndCirculatorsTags
*/
struct Forward_circulator_tag  : public virtual std::forward_circulator_tag{};

/*!
\ingroup PkgHandlesAndCirculatorsTags
*/
  struct Bidirectional_circulator_tag  : public virtual std::bidirectional_iterator_tag {};

/*!
\ingroup PkgHandlesAndCirculatorsTags
*/
  struct Random_access_circulator_tag : public virtual  std::random_access_circulator_tag {};




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
\ingroup PkgHandlesAndCirculatorsBaseClasses
*/
template <class T, class Dist, class Size>
struct Forward_circulator_base {};

/*!
\ingroup PkgHandlesAndCirculatorsBaseClasses
*/
template <class T, class Dist, class Size>
struct Bidirectional_circulator_base {};

/*!
\ingroup PkgHandlesAndCirculatorsBaseClasses
*/
template <class T, class Dist, class Size>
struct Random_access_circulator_base {};




/*!
\ingroup PkgHandlesAndCirculatorsBaseClasses
forward circulator.
*/
template <class T, class Dist, class Size>
class Forward_circulator_ptrbase {};

/*!
\ingroup PkgHandlesAndCirculatorsBaseClasses
bidirectional circulator.
*/
template <class T, class Dist, class Size>
class Bidirectional_circulator_ptrbase {};

/*!
\ingroup PkgHandlesAndCirculatorsBaseClasses
random access circulator.
*/
template <class T, class Dist, class Size>
class Random_access_circulator_ptrbase {};


} /* end namespace CGAL */
