namespace CGAL {

/*!
\ingroup PkgGenerators

computes a random convex `n`-gon by writing its vertices (oriented
counterclockwise) to `o`. The resulting polygon is scaled such
that it fits into the bounding box as specified by `pg`. Therefore
we cannot easily describe the resulting distribution.
\pre \f$ n \ge3\f$.

The function `random_convex_set_2` computes a random convex planar 
point set of given size where the points are drawn from a specific 
domain. 

Requirements 
-------------- 

<OL> 
<LI>`PointGenerator` is a model of the concept PointGenerator 
<LI>`Traits` is a model of the concept RandomConvexSetTraits_2 
<LI>`Point_generator::value_type` is equivalent to 
`Traits::Point_2` and `OutputIterator::value_type`. 
<LI>if `Traits` is not specified, 
`Point_generator::value_type` must be `Point_2< 
R >` for some representation class `R`, 
</OL> 

The default traits class `Default_traits` is 
`Random_convex_set_traits_2`. 
. 

\sa `CGAL::Random_points_in_square_2<Point_2, Creator>` 
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>` 

Implementation 
-------------- 

The implementation uses the centroid method 
described in \cite cgal:s-zkm-96 and has a worst case running time of \f$ O(r 
\cdot n + n \cdot \log n)\f$, where \f$ r\f$ is the time needed by `pg` 
to generate a random point. 

Example 
-------------- 

The following program displays a random convex 500-gon where the 
points are drawn uniformly from the unit square centered at the 
origin. 

\cgalexample{random_convex_set.cpp} 

*/
template < class OutputIterator, class PointGenerator,
class Traits > OutputIterator random_convex_set_2( std::size_t n,
OutputIterator o, const PointGenerator& pg, Traits t =
Default_traits);

} /* namespace CGAL */
