namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

The class `Delaunay_triangulation_on_sphere_2` is designed to represent
the Delaunay triangulation of a set of points on the 2-sphere.

A Delaunay triangulation of a set of points is a triangulation of the set of points that fulfills
the following <I>empty circle property</I> (also called <I>Delaunay property</I>): the circumscribing
sphere of any simplex of the triangulation contains no point of the set in its interior.
For a point set with no case of co-circularity of more than three points,
the Delaunay triangulation is unique, and defined as the dual of the Voronoi diagram of the points.

The setting of 3D points on the 2-sphere is particular in that the empty circle property
can be reduced to a single 3D orientation test \cgalCite{cgal:ccplr-redtp-10}.

\tparam Traits is the geometric traits; it must be a model of `DelaunayTriangulationOnSphereTraits_2`.

\tparam TDS is the triangulation data structure, which must be a model of `TriangulationDataStructure_2`
        whose vertex base must be a model of `TriangulationOnSphereVertexBase_2` and whose face base
        must be a model of `TriangulationOnSphereFaceBase_2`. \cgal provides a default instantiation
        for this parameter, which is the class
        `CGAL::Triangulation_data_structure_2 < CGAL::Triangulation_on_sphere_vertex_base_2<Traits>,
                                                CGAL::Triangulation_on_sphere_face_base_2<Traits> >`.

\sa `CGAL::Triangulation_on_sphere_2<Traits, TDS>`
*/
template< typename Traits, typename TDS >
class Delaunay_triangulation_on_sphere_2
  : public Triangulation_on_sphere_2<Traits, TDS>
{
public:
  /// \name Types
  /// @{

  /*!
  The traits class.
  */
  typedef Traits Geom_traits;

  /*!
  The point type representing a point on the sphere.
  */
  typedef Traits::Point_on_sphere_2 Point;

  /// @}

  /// \name Creation
  /// @{

  /*!
  introduces an empty triangulation and sets the center and radius of the sphere to `c` and `r` respectively.
  */
  Delaunay_triangulation_on_sphere_2(const Point_3& c, const FT r);

  /*!
  introduces an empty triangulation.

  \note The values for the center and radius must be either already set in the traits
  (if `gt` is provided) or must be set after the construction, using the function `set_center_and_radius()`.
  */
  Delaunay_triangulation_on_sphere_2(const Geom_traits& gt = Geom_traits());

  /*!
  introduces an empty triangulation, sets the center and radius of the sphere to `c` and `r` respectively,
  and inserts the point range `[first; beyond[`.

  \tparam PointOnSphereIterator must be a model of `InputIterator` with value type `Point_on_sphere_2` or `Point_3`.
  */
  template <typename PointOnSphereIterator>
  Delaunay_triangulation_on_sphere_2(PointOnSphereIterator first, PointOnSphereIterator beyond,
                                     const Point_3& center, const FT radius);

  /*!
  introduces an empty triangulation whose center and radius are set according to values within the traits
  and inserts the point range `[first;beyond[`.

  \warning It is the user's responsibility to ensure that the center and radius are set as intended in `gt`.

  \tparam PointOnSphereIterator must be a model of `InputIterator` with value type `Point_on_sphere_2` or `Point_3`.
  */
  template <typename PointOnSphereIterator>
  Delaunay_triangulation_on_sphere_2(PointOnSphereIterator first, PointOnSphereIterator beyond,
                                     const Geom_traits& gt = Geom_traits());

  /*!
  Copy constructor. All the vertices and faces are duplicated.
  */
  Delaunay_triangulation_on_sphere_2(const Delaunay_triangulation_on_sphere_2<Traits,TDS> &tr);

  /// @}

public:
  /// \name Predicates
  ///
  /// @{

  /*!
  returns the side of `p` with respect to the circle circumscribing the triangle associated with `f`.
  */
  Oriented_side side_of_oriented_circle(Face_handle f, const Point& p) const;

  /// @}

public:
  /// \name Insertion and Removal
  ///
  /// Methods for the insertion and removal of points on the sphere.
  /// In the degenerate setting of co-circular points, the Delaunay triangulation is known
  /// not to be uniquely defined. In this case, \cgal chooses a particular Delaunay triangulation
  /// using a symbolic perturbation scheme \cgalCite{cgal:dt-pvr3d-03}.
  ///
  /// @{

  /*!
  inserts the point `p`.
  If the point `p` coincides with an already existing vertex, this vertex is returned
  and the triangulation remains unchanged.
  The optional parameter `f` is used to give a hint about the location of `p`.
  */
  Vertex_handle insert(const Point& p, Face_handle f = Face_handle());

  /*!
  inserts the point `p` at the location given by `(lt, loc, li)`.

  \sa `Triangulation_on_sphere_2::locate()`
  */
  Vertex_handle insert(const Point& p, Locate_type& lt, Face_handle loc, int li );

  /*!
  Equivalent to `insert(p)`.
  */
  Vertex_handle push_back(const Point& p);

  /*!
  inserts the points in the range `[first, beyond)` and returns the number of inserted points.

  \tparam PointOnSphereIterator must be a model of `InputIterator` with value type `Point`.
  */
  template <class PointOnSphereIterator>
  size_type insert(PointOnSphereIterator first, PointOnSphereIterator beyond);

  /*!
  removes the vertex `v` from the triangulation.
  */
  void remove(Vertex_handle v);

  /// @}

  /// \name Queries
  /// @{

  /*!
  outputs the faces and boundary edges of the conflict zone of point `p` into output iterators.

  This function outputs in the container pointed to by `fit` the faces which are in conflict with point `p`,
  i. e., the faces whose circumcircle contains `p`. It outputs in the container pointed to by `eit`
  the boundary of the zone in conflict with `p`. The boundary edges of the conflict zone
  are output in counter-clockwise order and each edge is described through its incident face
  which is not in conflict with `p`. The function returns in a `std::pair` the resulting output iterators.

  \tparam OutItFaces is an output iterator with `Face_handle` as value type.
  \tparam OutItBoundaryEdges is an output iterator with `Edge` as value type.

  \warning This function makes uses of the member `tds_data` (see the concept `TriangulationDSFaceBase_2`)
           of the face to mark faces in conflict. It is the responsibility of the user to make sure the flags are cleared.

  \pre `dimension() == 2`.
  */
  template <typename OutputItFaces, typename OutputItBoundaryEdges>
  std::pair<OutputItFaces, OutputItBoundaryEdges>
  get_conflicts_and_boundary(const Point& p,
                             OutputItFaces fit,
                             OutputItBoundaryEdges eit,
                             Face_handle start) const;

  /// @}

  /// \name Voronoi Diagram
  ///
  /// The following member functions provide the elements of the dual Voronoi diagram.
  /// Two different embeddings are possible: a "straight" embedding, using line segments living
  /// in Euclidean 3D sphere, and a "curved" embedding, using arc segments on the sphere.
  ///
  /// Note that the following operations are constructions, which should be kept in mind in the choice
  /// of the underlying kernel.
  ///
  /// @{

  // Straight

  /*!
  returns the center of the circle circumscribed to face `f`
  */
  Point_3 dual(const Face_handle f) const;

  /*!
  returns the line segment with endpoints the circumcenters of the faces incident to the edge `e`.
  */
  Segment_3 dual(const Edge& e) const;

  /*!
  returns the line segment with endpoints the circumcenters of the faces incident to the edge `*ec`.
  */
  Segment_3 dual(const Edge_circulator ec) const;

  /*!
  returns the line segment with endpoints the circumcenters of the faces incident to the edge `*ei`.
  */
  Segment_3 dual(const All_edges_iterator ei) const;

  // Curved

  /*!
  returns the intersection of the dual of the face `f` and the sphere.

  \pre `dimension() == 2` and `f` is a solid face.
  */
  Point dual_on_sphere(const Face_handle f) const;

  /*!
  returns the arc of great circle with endpoints the circumcenters of the faces incident to the edge `e`.

  \pre `dimension() == 2` and `e` is not a ghost edge.
  */
  Arc_on_sphere_2 dual_on_sphere(const Edge& e) const;

  /*!
  returns the arc of great circle with endpoints the circumcenters of the faces incident to the edge `*ec`.

  \pre `dimension() == 2` and `*ec` is not a ghost edge.
  */
  Arc_on_sphere_2 dual_on_sphere(const Edge_circulator ec) const;

  /*!
  returns the arc of great circle with endpoints the circumcenters of the faces incident to the edge `*ei`.

  \pre `dimension() == 2` and `*ei` is not a ghost edge.
  */
  Arc_on_sphere_2 dual_on_sphere(const All_edges_iterator ei) const;

  /// @}

  /// \name Miscellaneous
  ///
  /// @{

  /*!
  tests the validity of the triangulation as a `Triangulation_on_sphere_2`
  and additionally tests the Delaunay property. This method is mainly useful
  for debugging Delaunay triangulation algorithms.
  */
  bool is_valid(bool verbose = false, int level = 0) const;

  /// @}

};

} // namespace CGAL
