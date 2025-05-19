namespace CGAL {

/*!
\ingroup PkgGeneratorsRef

\brief computes a random convex planar point set of given size where the
points are drawn from a specific domain.

The function computes a random convex `n`-gon by writing its vertices (oriented
counterclockwise) to `o`. The resulting polygon is scaled such
that it fits into the bounding box as specified by `pg`. Therefore
we cannot easily describe the resulting distribution.
\pre \f$ n \ge3\f$.



\cgalHeading{Requirements}

- `PointGenerator` is a model of the concept PointGenerator
- `Traits` is a model of the concept RandomConvexSetTraits_2
- `Point_generator::value_type` is equivalent to
`Traits::Point_2` and `OutputIterator::value_type`.
- if `Traits` is not specified,
`Point_generator::value_type` must be `Point_2<
R >` for some representation class `R`,


\sa `CGAL::Random_points_in_square_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>`

\cgalHeading{Implementation}

The implementation uses the centroid method
described in \cgalCite{cgal:s-zkm-96} and has a worst case running time of \cgalBigO{r
\cdot n + n \cdot \log n}, where \f$ r\f$ is the time needed by `pg`
to generate a random point.

\cgalHeading{Example}

The following program displays a random convex 500-gon where the
points are drawn uniformly from the unit square centered at the
origin.

\cgalExample{Generator/random_convex_set.cpp}

*/
template < class OutputIterator, class PointGenerator, class Traits >
OutputIterator random_convex_set_2( std::size_t n,
OutputIterator o, const PointGenerator& pg, Traits t = Random_convex_set_traits_2);

} /* namespace CGAL */
