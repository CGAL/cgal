namespace CGAL {

/*!
\ingroup PkgTriangulationOnSphere2TriangulationClasses

The class `Delaunay_triangulation_on_sphere_2` is designed to represent
the Delaunay triangulation of a set of points on the 2-sphere.

A Delaunay triangulation of a set of points is a triangulation of the sets of points that fulfills
the following <I>empty circle property</I> (also called <I>Delaunay property</I>): the circumscribing
sphere of any simplex of the triangulation contains no point of the set in its interior.
For a point set with no case of co-circularity of more than three points,
the Delaunay triangulation is unique, it is the dual of the Voronoi diagram of the points.

The setting of 3D points on the 2-sphere is particular in that the empty circle property
can be reduced to a single 3D orientation test \cgalCite{cgal:ccplr-redtp-10}.

\tparam Traits is the geometric traits, which must be a model of `DelaunayTriangulationTraits_2`.

\tparam TDS is the triangulation data structure, which must be a model of `TriangulationDataStructure_2`.
        \cgal provides a default instantiation for this parameter, which is the class
        `CGAL::Triangulation_data_structure_2 < CGAL::Triangulation_on_sphere_vertex_base_2<Traits>,
                                                CGAL::Triangulation_on_sphere_face_base_2<Traits> >`.

\sa `CGAL::Triangulation_on_sphere_2<Traits, TDS>`
*/
template< typename Traits, typename TDS >
class Delaunay_triangulation_on_sphere_2
  : public Triangulation_on_sphere_2<Traits, TDS>
{
public:
  /// \name Creation
  /// @{

  /*!
  Introduces an empty triangulation.
  */
  Delaunay_triangulation_on_sphere_2(const Traits& gt = Traits());

  /*!
  Introduces an empty triangulation and sets the center and radius of the sphere to `c` and `r` respectively.
  */
  Delaunay_triangulation_on_sphere_2(const Point_3& c, const FT r);

  /*!
  Copy constructor. All the vertices and faces are duplicated.
  */
  Delaunay_triangulation_on_sphere_2(const Delaunay_triangulation_on_sphere_2<Traits,TDS> &tr);

  /*!
  Equivalent to constructing an empty triangulation with the optional traits class argument
  and calling insert(first, last).
  */
  template <class InputIterator>
  Delaunay_triangulation_on_sphere_2(InputIterator first, InputIterator last, const Traits& gt = Traits());

  /// @}

public:
  /// \name Predicates
  ///
  /// @{

  /*!
  Returns the side of `p` with respect to the circle circumscribing the triangle associated with `f`.
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
  Inserts the point `p`.
  If the point `p` coincides with an already existing vertex, this vertex is returned
  and the triangulation remains unchanged.
  The optional parameter `f` is used to give a hint about the location of `p`.
  */
  Vertex_handle insert(const Point& p, Face_handle f = Face_handle());

  /*!
  Inserts the point `p` at the location given by `(lt, loc, li)`.
  \sa `Triangulation_on_sphere_2::locate()`
  */
  Vertex_handle insert(const Point& p, Locate_type& lt, Face_handle loc, int li );

  /*!
  Equivalent to `insert(p)`.
  */
  Vertex_handle push_back(const Point& p);

  /*!
  Inserts the points in the range `[first, last)` and returns the number of inserted points.
  \tparam PointInputIterator must be an input iterator with the value type `Point`.
  */
  template <class PointInputIterator>
  std::ptrdiff_t insert(PointInputIterator first, PointInputIterator last);

  /*!
  Removes the vertex `v` from the triangulation.
  */
  void remove(Vertex_handle v);

  /// @}

  /// \name Queries
  /// @{

  /*!
  Outputs the faces and boundary edges of the conflict zone of point `p` into output iterators.

  This function outputs in the container pointed to by `fit` the faces which are in conflict with point `p`,
  i. e., the faces whose circumcircle contains `p`. It outputs in the container pointed to by `eit` the
  the boundary of the zone in conflict with `p`. The boundary edges of the conflict zone
  are output in counter-clockwise order and each edge is described through its incident face
  which is not in conflict with `p`. The function returns in a `std::pair` the resulting output iterators.

  \tparam OutItFaces is an output iterator with `Face_handle` as value type.
  \tparam OutItBoundaryEdges is an output iterator with `Edge` as value type.

  \pre `dimension()==2`.
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
  /// Two different embeddings are possible: a "straight" embedding, using line segment living
  /// in Euclidean 3D sphere, and a "curved" embedding, using arc segments on the sphere.
  ///
  /// @{

  // Straight

  /*!
  Returns the center of the circle circumscribed to face `f`
  */
  Point_3 dual(const Face_handle f) const;

  /*!
  Returns the line segment with endpoints the circumcenters of the faces incident to the edge `e`.
  */
  Segment_3 dual(const Edge& e) const;

  /*!
  Same as above.
  */
  Segment_3 dual(const Edge_circulator ec) const;

  /*!
  Same as above.
  */
  Segment_3 dual(const All_edges_iterator ei) const;

  // Curved

  /*!
  Returns the intersection of the dual of the face `f` and the sphere
  */
  Point dual_on_sphere(const Face_handle f) const;

  /*!
  Returns the arc of great circle with endpoints the circumcenters of the faces incident to the edge `e`.
  */
  Arc_segment arc_dual(const Edge& e) const;

  /*!
  Same as above.
  */
  Arc_segment arc_dual(const Edge_circulator ec) const;

  /*!
  Same as above.
  */
  Arc_segment arc_dual(const All_edges_iterator ei) const;

  /// @}

  /// \name Miscellaneous
  ///
  /// @{

  /*!
  Tests the validity of the triangulation as a `Triangulation_on_sphere_2`
  and additionally tests the Delaunay property. This method is mainly useful
  for debugging Delaunay triangulation algorithms designed by the user.
  */
  bool is_valid(bool verbose = false, int level = 0) const;

  /// @}

};

} // namespace CGAL
