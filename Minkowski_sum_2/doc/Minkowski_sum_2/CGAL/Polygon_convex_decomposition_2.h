namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2

\anchor mink_refGreene_decomp 

The `Greene_convex_decomposition_2` class implements the approximation algorithm of 
Greene for the decomposition of an input polygon into convex 
sub-polygons \cite g-dpcp-83. This algorithm takes \f$ O(n \log n)\f$ 
time and \f$ O(n)\f$ space, where \f$ n\f$ is the size of the input polygon, 
and outputs a decomposition whose size is guaranteed to be no more 
than four times the size of the optimal decomposition. 

The `Polygon_2` type defined by the class is simply 
`Polygon_2<Kernel,Container>`. The `Container` parameter 
is by default `std::vector<typename Kernel::Point_2>`. 

\models ::PolygonConvexDecomposition_2 

\sa `CGAL::greene_approx_convex_partition_2` 

*/
template< typename Kernel, typename Container >
class Greene_convex_decomposition_2 {
public:

/// @}

}; /* end Greene_convex_decomposition_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2

\anchor mink_refHM_decomp 

The `Hertel_Mehlhorn_convex_decomposition_2` class implements the approximation algorithm of Hertel 
and Mehlhorn for decomposing a polygon into convex 
sub-polygons \cite hm-ftsp-83. This algorithm constructs a 
triangulation of the input polygon and proceeds by removing 
unnecessary triangulation edges. Given the triangulation, the 
algorithm requires \f$ O(n)\f$ time and space to construct a convex 
decomposition (where \f$ n\f$ is the size of the input polygon), whose 
size is guaranteed to be no more than four times the size of the 
optimal decomposition. 

The `Polygon_2` type defined by the class is simply 
`Polygon_2<Kernel,Container>`. The `Container` parameter 
is by default `std::vector<typename Kernel::Point_2>`. 

\models ::PolygonConvexDecomposition_2 

\sa `CGAL::approx_convex_partition_2` 

*/
template< typename Kernel, typename Container >
class Hertel_Mehlhorn_convex_decomposition_2 {
public:

/// @}

}; /* end Hertel_Mehlhorn_convex_decomposition_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2

\anchor mink_refopt_decomp 

The `Optimal_convex_decomposition_2` class provides an implementation of Greene's 
dynamic programming algorithm for optimal decomposition of a 
polygon into convex sub-polygons \cite g-dpcp-83. Note that 
this algorithm requires \f$ O(n^4)\f$ time and \f$ O(n^3)\f$ space in 
the worst case, where \f$ n\f$ is the size of the input polygon. 

The `Polygon_2` type defined by the class is simply 
`Polygon_2<Kernel,Container>`. The `Container` parameter 
is by default `std::vector<typename Kernel::Point_2>`. 

\models ::PolygonConvexDecomposition_2 

\sa `CGAL::optimal_convex_partition_2` 

*/
template< typename Kernel, typename Container >
class Optimal_convex_decomposition_2 {
public:

/// @}

}; /* end Optimal_convex_decomposition_2 */
} /* end namespace CGAL */
