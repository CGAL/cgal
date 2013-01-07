
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsTri3

The class `Kinetic::Regular_triangulation_3` maintains a triangulation of set of moving 
weighted points. Its interface is the same as 
`Kinetic::Delaunay_triangulation_3<Traits, Visitor, Triangulation>`. 

Note that the regular triangulation tracks as points are added to the `Kinetic::ActiveObjectsTable`, but not removed from it. 

The optional `Triangulation` template argument must be a model of 
`CGAL::RegularTriangulation_3` which has 
`Kinetic::Regular_triangulation_cell_base_3<Traits, Base>` as a 
cell base and 
`Kinetic::Regular_triangulation_vertex_base_3<Traits, Base>` as a 
vertex base. 

\sa `Kinetic::Delaunay_triangulation_3<Traits, Visitor, Triangulation>`
\sa `Kinetic::RegularTriangulationVisitor_3`

\cgalHeading{Example}

\cgalExample{Kinetic_data_structures/Kinetic_regular_triangulation_3.cpp} 

*/
template< typename Traits, typename Visitor, typename Triangulation >
class Regular_triangulation_3 {
public:

/// @}

}; /* end Kinetic::Regular_triangulation_3 */
} /* end namespace Kinetic */
} /* end namespace CGAL */
