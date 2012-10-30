/*!
\ingroup PkgMinkowskiSum2Concepts
\cgalconcept

A model of the `PolygonConvexDecomposition_2` concept is capable of decomposing an input 
polygon \f$ P\f$ into a set of convex sub-polygons \f$ P_1, \ldots, P_k\f$, such 
that \f$ \cup_{i=1}^{k}{P_k} = P\f$. 

\hasModel `CGAL::Small_side_angle_bisector_decomposition_2<Kernel,Container>`
\hasModel `CGAL::Optimal_convex_decomposition_2<Kernel,Container>`
\hasModel `CGAL::Hertel_Mehlhorn_convex_decomposition_2<Kernel,Container>` 
\hasModel `CGAL::Greene_convex_decomposition_2<Kernel,Container`

*/

class PolygonConvexDecomposition_2 {
public:

/// \name Types 
/// @{

/*! 
the geometric kernel type. 
*/ 
typedef Hidden_type Kernel; 

/*! 
the point type, used to represent polygon vertices. 
*/ 
typedef Hidden_type Point_2; 

/*! 
the polygon type. 
*/ 
typedef Hidden_type Polygon_2; 

/// @} 

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
PolygonConvexDecomposition_2 (); 

/// @} 

/// \name Operations 
/// @{

/*! 
decomposes the input polygon `P` into convex sub-polygons, 
and writes them to the output iterator `oi`. The value-type of the 
output iterator must be `Polygon_2`. 
The function returns a past-the-end iterator for the convex sub-polygons. 
*/ 
template <class OutputIterator> 
OutputIterator operator() (const Polygon_2& P, 
OutputIterator oi) const; 

/// @}

}; /* end PolygonConvexDecomposition_2 */
