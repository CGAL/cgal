 
namespace CGAL {

/*!
\ingroup PkgAdvancingFrontSurfaceReconstruction

The class `Advancing_front_surface_reconstruction` enables advanced users to provide the unstructured 
point cloud in a 3D Delaunay triangulation. The reconstruction algorithm then marks vertices and faces
in the triangulation as being on the 2D surface embedded in 3D space, and constructs a 2D triangulation
data structure that describes the surface.  The vertices and facets of the 2D triangulation data structure 
store handles to the vertices and faces of the 3D triangulation, which enables the user to explore the
2D as well as 3D neighborhood of vertices and facets of the surface. 

\tparam Gt must be a model of `Kernel`.
\tparam Dt must be a `Delaunay_triangulation_3` with
`Advancing_front_surface_reconstruction_vertex_base_3` and `Advancing_front_surface_reconstruction_cell_base_3` blended into the vertex and cell type, and the geometric traits class must be the `Exact_predicates_inexact_constructions_kernel`.

\tparam Filter must be a functor with `bool operator()(Gt::Point_3,Gt::Point_3,Gt::Point_3)` that allows the user to filter candidate triangles, for example based on its size.
        It defaults to a functor that always returns `false`.

*/

  template< typename Gt, typename Dt, typename Filter>
class Advancing_front_surface_reconstruction {
public:

/// \name Types 
/// @{

/*! 
  The type of the 2D triangulation data structure describing the reconstructed surface, being a model of `TriangulationDataStructure_2`.
  - The type `Triangulation_data_structure_2::Vertex` is model of the concept `TriangulationDataStructure_2::Vertex` and has additionally the 
  method `vertex_3()` that returns a `#Vertex_handle` to the associated 3D vertex.
  - The type `Triangulation_data_structure_2::Face` is model of the concept `TriangulationDataStructure_2::Face` and  has additionally the
  method `facet()` that returns the associated `#Facet`, and a method `bool is_on_surface()`
  for testing if a face is part of the reconstructed surface or a face incident to a boundary edge. 
 
  In case the surface has boundaries, the 2D surface has one vertex which is associated to the infinite
  vertex of the 3D triangulation.  
*/ 
  typedef unspecified_type Triangulation_data_structure_2; 


/*! 
The type of the 3D triangulation.

*/ 
  typedef Dt Triangulation_3; 

/*! 
The point type.

*/ 
  typedef typename Triangulation_3::Point Point;
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
  A bidirectional iterator range which enables to enumerate all points that were removed
  from the 3D Delaunay triangulation during the surface reconstruction. The value type
  of the iterator is `#Point`. 
*/ 
typedef unspecified_type Outlier_range; 

#if 0
/*! 
  A forward iterator which enables to visit all boundaries. It 
  visits the entry point of each boundary twice. This allows to
  detect that the traversal of a boundary is finished. One more increment
  brings us to the vertex on the next boundary. 
  The value type of the iterator is `#Vertex_handle`.
*/ 

#endif
  /*!
    A bidirectional iterator range which enables to visit all boundaries.
     The value type of the iterator is `Vertex_on_boundary_range`.
   */
typedef unspecified_type Boundary_range; 

 /*!
   A bidirectional iterator range which enables to visit all vertices on a boundary.
     The value type of the iterator is  `#Vertex_handle`
   */
typedef unspecified_type Vertex_on_boundary_range; 


/// @} 

/// \name Creation 
/// @{

/*! 
Constructor for the unstructured point cloud given as 3D Delaunay triangulation. 
*/ 
Advancing_front_surface_reconstruction(Dt& dt); 


/// @} 

/// \name Operations 
/// @{

/*! 
runs the surface reconstruction function.

\param radius_ratio_bound candidates incident to surface triangles which are not in the beta-wedge
       are discarded, if the ratio of their radius and the radius of the surface triangle is larger than `radius_ratio_bound`. 
       Described in Section \ref AFSR_Boundaries
\param beta half the angle of the wedge in which only the radius of triangles counts for the plausibility of candidates. 
       Described in Section \ref AFSR_Selection

*/ 
  void run(double radius_ratio_bound =5 , double beta = 0.52);

/*! 
returns the reconstructed surface.
*/ 
const Triangulation_data_structure_2& 
triangulation_data_structure_2() const; 

/*! 
returns the underlying 3D Delaunay triangulation. 
*/ 
const Triangulation_3& 
triangulation_3(); 




/*! 
returns an iterator range over the outliers.
*/ 
Outlier_range outliers(); 


/*! 
returns an iterator range over the boundaries. 
*/ 
Boundary_range boundaries(); 


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
returns `true` if the facet `f2` is on the surface.
*/ 
bool
has_on_surface(Triangulation_data_structure_2::Face_handle f2) const; 

/*! 
returns `true` if the vertex `v2` is on the surface.
*/ 
bool
has_on_surface(Triangulation_data_structure_2::Vertex_handle v2) const;
/// @} 




}; /* end Advancing_front_surface_reconstruction */


/*!
\ingroup PkgAdvancingFrontSurfaceReconstruction

For a sequence of points computes a sequence of index triples
describing the faces of the reconstructed surface.

\tparam PointInputIterator must be an input iterator with 3D points from the `Exact_predicates_inexact_constructions_kernel` as value type.
\tparam IndicesOutputIterator must be an output iterator to which 
`CGAL::cpp11::tuple<std::size_t,std::size_t,std::size_t>` can be assigned.

\param b iterator on the first point of the sequence
\param e past the end iterator of the point sequence
\param out output iterator
\param radius_ratio_bound candidates incident to surface triangles which are not in the beta-wedge
       are discarded, if the ratio of their radius and the radius of the surface triangle is larger than `radius_ratio_bound`. 
       Described in Section \ref AFSR_Boundaries
\param beta half the angle of the wedge in which only the radius of triangles counts for the plausibility of candidates. 
       Described in Section \ref AFSR_Selection

*/
  template <class PointInputIterator, IndicesOutputIterator>
  IndicesOutputIterator advancing_front_surface_reconstruction(PointInputIterator b, PointInputIterator e, IndicesOutputIterator out, double radius_ratio_bound = 5, double beta= 0.52 );


/*!
\ingroup PkgAdvancingFrontSurfaceReconstruction

For a sequence of points computes a sequence of index triples
describing the faces of the reconstructed surface.

\tparam PointInputIterator must be an input iterator with 3D points from the `Exact_predicates_inexact_constructions_kernel` as value type.
\tparam IndicesOutputIterator must be an output iterator to which 
`CGAL::cpp11::tuple<std::size_t,std::size_t,std::size_t>` can be assigned.
\tparam Filter must be a functor with `bool operator()(Gt::Point_3,Gt::Point_3,Gt::Point_3)`.

\param b iterator on the first point of the sequence
\param e past the end iterator of the point sequence
\param out output iterator
\param radius_ratio_bound candidates incident to surface triangles which are not in the beta-wedge
       are discarded, if the ratio of their radius and the radius of the surface triangle is larger than `radius_ratio_bound`. 
       Described in Section \ref AFSR_Boundaries
\param beta half the angle of the wedge in which only the radius of triangles counts for the plausibility of candidates. 
       Described in Section \ref AFSR_Selection
\param filter allows the user to filter candidate triangles, for example based on their size.   

*/
  template <class PointInputIterator, IndicesOutputIterator, class Filter>
  IndicesOutputIterator advancing_front_surface_reconstruction(PointInputIterator b, PointInputIterator e, IndicesOutputIterator out, Filter filter, double radius_ratio_bound = 5, double beta= 0.52 );



} /* end namespace CGAL */
