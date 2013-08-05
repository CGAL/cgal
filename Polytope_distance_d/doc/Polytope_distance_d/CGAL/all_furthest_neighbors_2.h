namespace CGAL {

/*!
\ingroup PkgOptimalDistances

\brief computes all furthest neighbors for the vertices of the convex
polygon described by the range [`points_begin`, `points_end`), writes
their indices (relative to `points_begin`) to `o`\cgalFootnote{the
furthest neighbor of `points_begin[i]` is `points_begin[i-th number
written to o]`} and returns the past-the-end iterator of this
sequence.

The function `all_furthest_neighbors_2()` computes all furthest 
neighbors for the vertices of a convex polygon \f$ P\f$, i.e.\ for each 
vertex \f$ v\f$ of \f$ P\f$ a vertex \f$ f_v\f$ of \f$ P\f$ such that the distance 
between \f$ v\f$ and \f$ f_v\f$ is maximized. 

\pre The points denoted by the non-empty range [`points_begin`, `points_end`) form the boundary of a convex polygon \f$ P\f$ (oriented clock- or counterclockwise). 

The geometric types and operations to be used for the computation 
are specified by the traits class parameter `t`. This parameter 
can be omitted if `RandomAccessIterator` refers to a point type 
from a `Kernel`. In this case, the kernel is used as default 
traits class. 

\cgalHeading{Requires}

<OL> 
<LI>If `t` is specified explicitly, `Traits` is a model 
for `AllFurthestNeighborsTraits_2`. 
<LI>Value type of `RandomAccessIterator` is 
`Traits::Point_2` or, if `t` is not specified explicitly,
`K::Point_2` where `K` is a model for `Kernel`. 
<LI>`OutputIterator` accepts `int` as value type. 
</OL> 

\sa `AllFurthestNeighborsTraits_2` 
\sa `CGAL::monotone_matrix_search()` 

\cgalHeading{Implementation}

The implementation uses monotone matrix search \cgalCite{akmsw-gamsa-87}. 
Its runtime complexity is linear in the number of
vertices of \f$ P\f$.

\cgalHeading{Example}

The following code generates a random convex polygon 
`p` with ten vertices, computes all furthest neighbors and 
writes the sequence of their indices (relative to 
`points_begin`) to `cout` (e.g. a sequence of 
`4788911224` means the furthest neighbor of 
`points_begin[0]` is `points_begin[4]`, the furthest 
neighbor of `points_begin[1]` is `points_begin[7]` etc.). 

\cgalExample{Polytope_distance_d/all_furthest_neighbors_2.cpp} 

*/
template < class RandomAccessIterator, class
OutputIterator, class Traits > OutputIterator
all_furthest_neighbors_2( RandomAccessIterator points_begin,
RandomAccessIterator points_end, OutputIterator o, Traits t =
Default_traits);

} /* namespace CGAL */

