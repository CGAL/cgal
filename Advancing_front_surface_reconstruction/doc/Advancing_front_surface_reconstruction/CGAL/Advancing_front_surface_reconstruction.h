 
namespace CGAL {

/*!
\ingroup PkgAdvancingFrontSurfaceReconstruction

The class `Advancing_front_surface_reconstruction` 

\tparam Dt must be a `Delaunay_triangulation_3` with
`Advancing_front_surface_reconstruction_vertex_base_3` and `Advancing_front_surface_reconstruction_cell_base_3` blended into the vertex and cell type, and the geometric traits class must be the `Exact_predicates_inexact_constructions_kernel`.

*/

template< typename Dt>
class Advancing_front_surface_reconstruction {
public:

/// \name Types 
/// @{

/*! 
  The type of the 2D triangulation data structure describing the reconstructed surface.
  The type `TDS_2::Vertex` is model of the concept `TriangulationDataStructure_2::Vertex` and has additionally the 
  method `vertex_3()` that returns a `Vertex_handle` to the associated 3D vertex
  The type `TDS_2::Face` is model of the concept `TriangulationDataStructure_2::Face` and  has additionally the
  method `facet()` that returns the associated `Facet`, and a method `bool is_on_surface()`
  for testing if a face is part of the reconstructed surface or a face incident to a boundary edge.  
  In case the surface has boundaries, the 2D surface has one vertex which is associated to the infinite
  vertex of the 3D triangulation.  
*/ 
  typedef unspecified_type TDS_2; 

/*! 
The type of the 3D triangulation.

*/ 
  typedef Dt Triangulation_3; 

/*! 
The vertex handle type of the 3D triangulation.

*/ 
  typedef typename Triangulation_3::Vertex_handle Vertex_handle;

/*! 
The cell handle type of the 3D triangulation.

*/ 
  typedef typename Triangulation_3::Cell_handle Cell_handle;

/*! 
The facet type of the 3D triangulation.

*/ 
  typedef typename Triangulation_3::Facet Facet; 


/*! 
  A bidirectional iterator which allows to enumerate all points that were removed
  from the 3D Delaunay triangulation during the surface reconstruction. The value type
  of the iterator is `Triangulation_3::Point_3`. 
*/ 
typedef unspecified_type Outlier_iterator; 

/*! 
  A forward iterator which allows to visit all boundaries. It 
  visits the entry point of each boundary twice. This allows to
  detect that the traversal of a boundary is finished. One more increment
  brings us to the vertex on the next boundary. 
  The value type of the iterator is `Triangulation_3::Vertex_handle`.
*/ 
typedef unspecified_type Boundary_iterator; 


/// @} 

/// \name Creation 
/// @{

/*! 
Initializes from a 3D Delaunay triangulation of a point set. 
*/ 
Advancing_front_surface_reconstruction(Dt& dt); 


/// @} 

/// \name Operations 
/// @{

/*! 
calls the surface reconstruction function with the default parameters.

\param K  described in Section \ref AFSR_Boundaries
\param beta described in Section \ref AFSR_Selection

*/ 
  void operator()(double K=5, double beta= 0.18);

/*! 
returns the reconstructed surface.
*/ 
const TDS_2& 
tds_2(); 

/*! 
returns the underlying 3D Delaunay triangulation. 
*/ 
const Triangulation_3& 
triangulation_3(); 




/*! 
An iterator over the outliers.
*/ 
Outlier_iterator outliers_begin(); 

/*! 
Past-the-end iterator. 
*/ 
Outlier_iterator outliers_end(); 

/*! 
An iterator over the boundary vertices. 
*/ 
Boundary_iterator boundaries_begin(); 

/*! 
Past-the-end iterator. 
*/ 
Boundary_iterator boundaries_end(); 

/// @} 

/// \name Predicates 
/// @{

/*! 
returns `true` if the reconstructed surface has boundaries. 
*/ 
bool 
has_boundaries() const; 


/*! 
returns `true` if the facet is on the surface.
*/ 
bool
has_on_surface(Facet f) const; 

/*! 
returns `true` if the facet f is on the surface.
*/ 
bool
has_on_surface(TDS_2::Face_handle f2) const; 

/*! 
returns `true` if the facet f is on the surface.
*/ 
bool
has_on_surface(TDS_2::Vertex_handle v2) const;
/// @} 




}; /* end Advancing_front_surface_reconstruction */


/*!
\ingroup PkgAdvancingFrontSurfaceReconstruction

For a sequence of points computes a sequence of index triples
describing the faces of the reconstructed surface.

\tparam PointInputIterator must be an input iterator with 3D points from the `Exact_predicates_inexact_constructions_kernel` as value type.
\tparam IndicesOutputIterator must be an output iterator to which 
`CGAL::cpp11::tuple<std::size_t,std::size_t,std::size_t>` can be assigned.

\param K  described in Section \ref AFSR_Boundaries
\param beta described in Section \ref AFSR_Selection

*/
  template <class PointInputIterator, IndicesOutputIterator>
  void advancing_front_surface_reconstruction(PointInputIterator b, PointInputIterator e, IndicesOutputIterator out, double K=5, double beta= 0.18 );


} /* end namespace CGAL */
