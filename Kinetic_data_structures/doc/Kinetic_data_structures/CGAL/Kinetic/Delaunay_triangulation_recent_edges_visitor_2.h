namespace CGAL { namespace Kinetic {

/*!
\ingroup PkgKdsTri2

The class `Kinetic::Delaunay_triangulation_recent_edges_visitor_2` provides a model of 
`Kinetic::DelaunayTriangulationVisitor_2` which tracks which edges were created in 
the most recent change. 

\cgalModels `Kinetic::DelaunayTriangulationVisitor_2`

\sa `Kinetic::Delaunay_triangulation_2<Traits, Triangulation, Visitor>` 

*/
template< typename Triangulation >
class Delaunay_triangulation_recent_edges_visitor_2 {
public:

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Delaunay_triangulation_recent_edges_visitor_2(); 

/*!
The iterator through the recently created edges. 
*/ 
typedef unspecified_type iterator; 

/// @} 

/// \name Operations 
/// @{

/*!
Begin iteration through the recent edges. 
*/ 
iterator begin() const; 

/*!
End iteration through the recent edges. 
*/ 
iterator end() const; 

/*!
Returns true if this edge exists in the set. 
*/ 
bool contains(Triangulation::Edge) const; 

/// @}

}; /* end Kinetic::Delaunay_triangulation_recent_edges_visitor_2 */
} /* end namespace Kinetic */
} /* end namespace CGAL */

