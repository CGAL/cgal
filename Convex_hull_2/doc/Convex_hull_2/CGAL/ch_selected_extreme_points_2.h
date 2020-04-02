namespace CGAL {

/*!
\ingroup PkgConvexHull2Extreme

The function `ch_e_point()` finds a point of a given set
of input points with maximal \f$ x\f$ coordinate.

It traverses the range [`first`,`beyond`).
After execution, the value of
`e` is an iterator in the range such that `*e` \f$ \ge_{xy}\f$
`*it` for all iterators `it` in the range.

The default traits class `Default_traits` is the kernel in which the
value type of `ForwardIterator` is defined.

\cgalHeading{Requirements}

`Traits` defines a type `Traits::Less_xy_2` as described in
the concept `ConvexHullTraits_2` and the corresponding member
function that returns an instance of this type.

\sa `CGAL::ch_nswe_point()`
\sa `CGAL::ch_n_point()`
\sa `CGAL::ch_ns_point()`
\sa `CGAL::ch_s_point()`
\sa `CGAL::ch_w_point()`
\sa `CGAL::ch_we_point()`


*/
template <class ForwardIterator>
void
ch_e_point( ForwardIterator first, ForwardIterator beyond,
ForwardIterator& e,
const Traits & ch_traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgConvexHull2Extreme

The function `ch_n_point()` finds a point in a given set
of input points with maximal \f$ y\f$ coordinate.

It traverses the range [`first`,`beyond`).
After execution, the value of
`n` is an iterator in the range such that `*n` \f$ \ge_{yx}\f$
`*it` for all iterators `it` in the range.

The default traits class `Default_traits` is the kernel in which the type
`ForwardIterator::value_type` is defined.

\cgalHeading{Requirements}

`Traits` defines the type `Traits::Less_yx_2` as specified in
the concept `ConvexHullTraits_2` and the corresponding member
function that returns an instance of this type.

\sa `CGAL::ch_e_point()`
\sa `CGAL::ch_nswe_point()`
\sa `CGAL::ch_ns_point()`
\sa `CGAL::ch_s_point()`
\sa `CGAL::ch_w_point()`
\sa `CGAL::ch_we_point()`


*/
template <class ForwardIterator>
void
ch_n_point( ForwardIterator first, ForwardIterator beyond,
ForwardIterator& n,
const Traits & ch_traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgConvexHull2Extreme

The function `ch_ns_point()` finds the points of a given set
of input points with minimal and maximal \f$ x\f$ coordinates.

It traverses the range [`first`,`beyond`).
After execution, the value of
`n` is an iterator in the range such that `*n` \f$ \ge_{yx}\f$
`*it` for all iterators `it` in the range. Similarly, for
`s` the inequality `*s` \f$ \le_{yx}\f$ `*it`
holds for all iterators in the range.

The default traits class `Default_traits` is the kernel in which the
value type of `ForwardIterator` is defined.

\cgalHeading{Requirements}

`Traits` defines the type `Traits::Less_yx_2` as specified in
the concept `ConvexHullTraits_2` and the corresponding member
function that returns an instance of this type.

\sa `CGAL::ch_e_point()`
\sa `CGAL::ch_nswe_point()`
\sa `CGAL::ch_n_point()`
\sa `CGAL::ch_s_point()`
\sa `CGAL::ch_w_point()`
\sa `CGAL::ch_we_point()`


*/
template <class ForwardIterator>
void
ch_ns_point( ForwardIterator first, ForwardIterator beyond,
ForwardIterator& n,
ForwardIterator& s,
const Traits & ch_traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgConvexHull2Extreme

The function `ch_nswe_point()` finds the four extreme points of a given set
of input points using a linear scan of the input points.
That is, it determines the points with maximal \f$ y\f$, minimal \f$ y\f$,
minimal \f$ x\f$, and maximal \f$ x\f$ coordinates.

It traverses the range [`first`,`beyond`).
After execution, the value of
`n` is an iterator in the range such that `*n` \f$ \ge_{yx}\f$
`*it` for all iterators `it` in the range. Similarly, for
`s`, `w`, and `e` the inequalities `*s` \f$ \le_{yx}\f$
`*it`, `*w` \f$ \le_{xy}\f$ `*it`, and `*e`
\f$ \ge_{xy}\f$ `*it` hold for all iterators
`it` in the range.

\cgalHeading{Requirements}

`Traits` contains the following subset of types from
the concept `ConvexHullTraits_2` and their corresponding member
functions that return instances of these types:
<UL>
<LI>`Traits::Less_xy_2`,
<LI>`Traits::Less_yx_2`.
</UL>

The default traits class `Default_traits` is the kernel in which the
type `ForwardIterator::value_type` is defined.

\sa `CGAL::ch_e_point()`
\sa `CGAL::ch_n_point()`
\sa `CGAL::ch_ns_point()`
\sa `CGAL::ch_s_point()`
\sa `CGAL::ch_w_point()`
\sa `CGAL::ch_we_point()`


*/
template <class ForwardIterator>
void
ch_nswe_point( ForwardIterator first, ForwardIterator beyond,
ForwardIterator& n,
ForwardIterator& s,
ForwardIterator& w,
ForwardIterator& e,
const Traits & ch_traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgConvexHull2Extreme

The function `ch_s_point()` finds a points in a given set
of input points with minimal \f$ y\f$ coordinates.

It traverses the range [`first`,`beyond`).
After execution, the value of
`s` is an iterator in the range such that `*s` \f$ \le_{yx}\f$
`*it` for all iterators `it` in the range.

The default traits class `Default_traits` is the kernel in which the
type `ForwardIterator::value_type` is defined.

\cgalHeading{Requirements}

`Traits` defines the type `Traits::Less_yx_2` as specified in
the concept `ConvexHullTraits_2` and the corresponding member
function that returns an instance of this type.

\sa `CGAL::ch_e_point()`
\sa `CGAL::ch_nswe_point()`
\sa `CGAL::ch_n_point()`
\sa `CGAL::ch_ns_point()`
\sa `CGAL::ch_w_point()`
\sa `CGAL::ch_we_point()`


*/
template <class ForwardIterator>
void
ch_s_point( ForwardIterator first, ForwardIterator beyond,
ForwardIterator& s,
const Traits & ch_traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgConvexHull2Extreme

The function `ch_we_point()` finds two points of a given set
of input points with minimal and maximal \f$ x\f$ coordinates.

It traverses the range [`first`,`beyond`).
After execution, the value of
`w` is an iterator in the range such that `*w` \f$ \le_{xy}\f$
`*it` for all iterators `it` in the range. Similarly, for
`e` the inequality `*e` \f$ \ge_{xy}\f$ `*it`
holds for all iterators in the range.

The default traits class `Default_traits` is the kernel in which the
type `ForwardIterator::value_type` is defined.

\cgalHeading{Requirements}

`Traits` defines the type `Traits::Less_xy_2` as specified in
the concept `ConvexHullTraits_2` and the corresponding member
function that returns an instance of this type.

\sa `CGAL::ch_e_point()`
\sa `CGAL::ch_nswe_point()`
\sa `CGAL::ch_n_point()`
\sa `CGAL::ch_ns_point()`
\sa `CGAL::ch_s_point()`
\sa `CGAL::ch_w_point()`


*/
template <class ForwardIterator>
void
ch_we_point( ForwardIterator first, ForwardIterator beyond,
ForwardIterator& w,
ForwardIterator& e,
const Traits & ch_traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgConvexHull2Extreme

The function `ch_w_point()` finds a point in a given set
of input points with minimal \f$ x\f$ coordinate.

It traverses the range [`first`,`beyond`).
After execution, the value of
`w` is an iterator in the range such that `*w` \f$ \le_{xy}\f$
`*it` for all iterators `it` in the range.

\cgalHeading{Requirements}

`Traits` defines the type `Traits::Less_xy_2` as specified in
the concept `ConvexHullTraits_2` and the corresponding member
function that returns an instance of this type.

The default traits class `Default_traits` is the kernel in which the
type `ForwardIterator::value_type` is defined.

\sa `CGAL::ch_e_point()`
\sa `CGAL::ch_nswe_point()`
\sa `CGAL::ch_n_point()`
\sa `CGAL::ch_ns_point()`
\sa `CGAL::ch_s_point()`
\sa `CGAL::ch_we_point()`


*/
template <class ForwardIterator>
void
ch_w_point( ForwardIterator first, ForwardIterator beyond,
ForwardIterator& w,
const Traits & ch_traits = Default_traits);

} /* namespace CGAL */
