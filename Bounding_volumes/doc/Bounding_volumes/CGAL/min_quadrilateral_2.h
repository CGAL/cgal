namespace CGAL {

/*!
\ingroup PkgBoundingVolumes

The function computes a minimum area enclosing 
parallelogram \f$ A(P)\f$ of a given convex point set \f$ P\f$. Note that 
\f$ R(P)\f$ is not necessarily axis-parallel, and it is in general not 
unique. The focus on convex sets is no restriction, since any 
parallelogram enclosing \f$ P\f$ - as a convex set - contains the convex 
hull of \f$ P\f$. For general point sets one has to compute the convex hull 
as a preprocessing step. 

computes a minimum area enclosing parallelogram of the point set 
described by [`points_begin`, `points_end`), writes its 
vertices (counterclockwise) to `o` and returns the past-the-end 
iterator of this sequence. 
If the input range is empty, `o` remains unchanged. 

If the input range consists of one element only, this point is written 
to `o` four times. 

\pre The points denoted by the range [`points_begin`, `points_end`) form the boundary of a simple convex polygon \f$ P\f$ in counterclockwise orientation. 

The geometric types and operations to be used for the computation 
are specified by the traits class parameter `t`. The parameter can be 
omitted, if `ForwardIterator` refers to a two-dimensional point 
type from one the \cgal kernels. In this case, a default traits class 
(`Min_quadrilateral_default_traits_2<K>`) is used. 

<OL>
<LI>If `Traits` is specified, it must be a model for 
`MinQuadrilateralTraits_2` and the value type `VT` of 
`ForwardIterator` is `Traits::Point_2`. Otherwise 
`VT` must be `CGAL::Point_2<K>` for some kernel 
`K`. 
<LI>`OutputIterator` must accept `VT` as value type. 
</OL> 

\sa `CGAL::min_rectangle_2()` 
\sa `CGAL::min_strip_2()` 
\sa `MinQuadrilateralTraits_2` 
\sa `CGAL::Min_quadrilateral_default_traits_2<K>` 

\cgalHeading{Implementation}

We use a rotating caliper 
algorithm 
\cgalCite{stvwe-mepa-95}, \cgalCite{v-fmep-90} with worst case running time linear 
in the number of input points. 

\cgalHeading{Example}

The following code generates a random convex polygon 
`P` with 20 vertices and computes the minimum enclosing 
parallelogram of `P`. 

\cgalExample{Min_quadrilateral_2/minimum_enclosing_parallelogram_2.cpp} 

*/

template < class ForwardIterator, class OutputIterator, class Traits >
OutputIterator
min_parallelogram_2(
ForwardIterator points_begin,
ForwardIterator points_end,
OutputIterator o,
Traits& t = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgBoundingVolumes

The function computes a minimum area enclosing rectangle 
\f$ R(P)\f$ of a given convex point set \f$ P\f$. Note that \f$ R(P)\f$ is not 
necessarily axis-parallel, and it is in general not unique. The focus 
on convex sets is no restriction, since any rectangle enclosing 
\f$ P\f$ - as a convex set - contains the convex hull of \f$ P\f$. For general 
point sets one has to compute the convex hull as a preprocessing step. 

computes a minimum area enclosing rectangle of the point set described 
by [`points_begin`, `points_end`), writes its vertices 
(counterclockwise) to `o`, and returns the past-the-end iterator 
of this sequence. 

If the input range is empty, `o` remains unchanged. 

If the input range consists of one element only, this point is written 
to `o` four times. 

\pre The points denoted by the range [`points_begin`, `points_end`) form the boundary of a simple convex polygon \f$ P\f$ in counterclockwise orientation. 

The geometric types and operations to be used for the computation are 
specified by the traits class parameter `t`. The parameter can be 
omitted, if `ForwardIterator` refers to a two-dimensional point 
type from one the \cgal kernels. In this case, a default traits class 
(`Min_quadrilateral_default_traits_2<K>`) is used. 

<OL>
<LI>If `Traits` is specified, it must be a model for 
`MinQuadrilateralTraits_2` and the value type `VT` of 
`ForwardIterator` is `Traits::Point_2`. Otherwise `VT` 
must be `CGAL::Point_2<K>` for some kernel `K`. 
<LI>`OutputIterator` must accept `VT` as value type. 
</OL> 

\sa `CGAL::min_parallelogram_2()` 
\sa `CGAL::min_strip_2()` 
\sa `MinQuadrilateralTraits_2` 
\sa `CGAL::Min_quadrilateral_default_traits_2<K>` 

\cgalHeading{Implementation}

We use a rotating caliper 
algorithm \cgalCite{t-sgprc-83} 
with worst case running time linear in the number of input points. 

\cgalHeading{Example}

The following code generates a random convex polygon 
`P` with 20 vertices and computes the minimum enclosing 
rectangle of `P`. 

\cgalExample{Min_quadrilateral_2/minimum_enclosing_rectangle_2.cpp} 

*/

template < class ForwardIterator, class OutputIterator, class Traits >
OutputIterator
min_rectangle_2(
ForwardIterator points_begin,
ForwardIterator points_end,
OutputIterator o,
Traits& t = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgBoundingVolumes

The function computes a minimum width enclosing strip 
\f$ S(P)\f$ of a given convex point set \f$ P\f$. A strip is the closed region 
bounded by two parallel lines in the plane. Note that \f$ S(P)\f$ is not 
unique in general. The focus on convex sets is no restriction, since any 
parallelogram enclosing \f$ P\f$ - as a convex set - contains the convex 
hull of \f$ P\f$. For general point sets one has to compute the convex hull 
as a preprocessing step. 

computes a minimum enclosing strip of the point set described by 
[`points_begin`, `points_end`), writes its two bounding lines 
to `o` and returns the past-the-end iterator of this sequence. 

If the input range is empty or consists of one element only, `o` 
remains unchanged. 

\pre The points denoted by the range [`points_begin`, `points_end`) form the boundary of a simple convex polygon \f$ P\f$ in counterclockwise orientation. 

The geometric types and operations to be used for the computation 
are specified by the traits class parameter `t`. The parameter 
can be omitted, if `ForwardIterator` refers to a two-dimensional 
point type from one the \cgal kernels. In this case, a default 
traits class (`Min_quadrilateral_default_traits_2<K>`) is 
used. 

<OL>
<LI>If `Traits` is specified, it must be a model for 
`MinQuadrilateralTraits_2` and the value type `VT` of 
`ForwardIterator` is `Traits::Point_2`. Otherwise `VT` 
must be `CGAL::Point_2<K>` for some kernel `K`. 
<LI>`OutputIterator` must accept `Traits::Line_2` as value type. 
</OL> 

\sa `CGAL::min_rectangle_2()` 
\sa `CGAL::min_parallelogram_2()` 
\sa `MinQuadrilateralTraits_2` 
\sa `CGAL::Min_quadrilateral_default_traits_2<K>` 

\cgalHeading{Implementation}

We use a rotating caliper 
algorithm \cgalCite{t-sgprc-83} 
with worst case running time linear in the number of input points. 

\cgalHeading{Example}

The following code generates a random convex polygon 
`P` with 20 vertices and computes the minimum enclosing 
strip of `P`. 

\cgalExample{Min_quadrilateral_2/minimum_enclosing_strip_2.cpp} 

*/

template < class ForwardIterator, class OutputIterator, class Traits >
OutputIterator
min_strip_2(
ForwardIterator points_begin,
ForwardIterator points_end,
OutputIterator o,
Traits& t = Default_traits);

} /* namespace CGAL */

