// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

namespace CGAL
{

/*!
\ingroup PkgPeriodic2Triangulation2MainClasses

The class `Periodic_2_Delaunay_triangulation_2` represents a
Delaunay triangulation in two-dimensional periodic space.

\cgalHeading{Parameters}

The class `Periodic_2_Delaunay_triangulation_2` has two template parameters. The first one
\tparam Traits is the geometric traits, it is to be instantiated by a
model of the concept `Periodic_2DelaunayTriangulationTraits_2`.

\tparam Tds is the triangulation data data structure and must be a model of `TriangulationDataStructure_2`
whose vertex and face are models of `Periodic_2TriangulationVertexBase_2` and `Periodic_2TriangulationFaceBase_2`.
It defaults to:
\code
CGAL::Triangulation_data_structure_2<
  CGAL::Periodic_2_triangulation_vertex_base_2<Gt>,
  CGAL::Periodic_2_triangulation_face_base_2<Gt> > >
\endcode

\cgalHeading{Implementation}

Insertion is implemented by inserting in the triangulation, then
performing a sequence of Delaunay flips. The number of flips is \cgalBigO{d}
if the new vertex is of degree \f$ d\f$ in the new triangulation. For
points distributed uniformly at random, insertion takes time \cgalBigO{1} on
average.

Removal calls the removal in the triangulation and then
re-triangulates the hole in such a way that the Delaunay criterion is
satisfied. Removal of a vertex of degree \f$ d\f$ takes time \cgalBigO{d^2}. The
expected degree \f$ d\f$ is \cgalBigO{1} for a random vertex in the
triangulation.

After a point location step, the nearest neighbor is found in time
\cgalBigO{n} in the worst case, but in expected time \cgalBigO{1} on average for
vertices distributed uniformly at random and any query point.

\sa `CGAL::Periodic_2_triangulation_2<Traits,Tds>`
\sa `CGAL::Periodic_2_triangulation_hierarchy_2<Tr>`
\sa `CGAL::Delaunay_triangulation_2<Traits,Tds>`
*/
template< typename Traits, typename Tds >
class Periodic_2_Delaunay_triangulation_2 : public Periodic_2_triangulation_2<Traits, Tds>
{
public:

/// \name Creation
/// @{

  /*!
  Creates an empty periodic Delaunay triangulation `dt`, with
  `domain` as original domain and possibly specifying
  a traits class `traits`.
  \pre `domain` is a square.
  */
  Periodic_2_Delaunay_triangulation_2(
    const Iso_rectangle & domain = Iso_rectangle(0, 0, 1, 1),
    const Geom_traits & traits = Geom_traits());

  /*!
  Copy constructor.
  */
  Periodic_2_Delaunay_triangulation_2 (const
                                       Periodic_2_Delaunay_triangulation_2 & dt1);

  /*!
  Equivalent to constructing an empty triangulation with the optional
  domain and traits class arguments and calling `insert(first,last)`.
  \pre The `value_type` of `first` and `last` are `Point`s lying inside `domain` and `domain` is a square.
  */
  template < class InputIterator >
  Periodic_2_Delaunay_triangulation_2 (
    InputIterator first,
    InputIterator last,
    const Iso_rectangle & domain = Iso_rectangle(0, 0, 1, 1),
    const Geom_traits & traits = Geom_traits());

/// @}

/// \name Insertion and removal
/// The following methods insert points in the triangulation ensuring
/// the empty circle property of Delaunay triangulations. The inserted
/// points need to lie in the original domain (see Section \ref
/// P2Triangulation2secspace of the user manual). In the degenerate
/// case when there are co-circular points, the Delaunay triangulation
/// is known not to be uniquely defined. In this case, \cgal chooses a
/// particular Delaunay triangulation using a symbolic perturbation
/// scheme \cgalCite{cgal:dt-pvr3d-03}. Note that insertion of a new point
/// can cause a switch from computing in the 9-sheeted covering space
/// to computing in the 1-sheeted covering space, which invalidates
/// some `Vertex_handle`s and `Face_handle`s.
/// @{

  /*!
  Inserts point `p` in the triangulation and returns the
  corresponding vertex. The optional argument `start` is used as a
  starting place for the point location. \pre `p` lies in the original domain `domain`.
  */
  Vertex_handle insert(const Point & p, Face_handle start =
                         Face_handle() );

  /*!
  Inserts point `p` in the triangulation and returns the
  corresponding vertex. Similar to the above `insert()` function,
  but takes as additional parameter the return values of a previous
  location query. See description of
  `Periodic_2_triangulation_2::locate()`. \pre `p` lies in the original domain `domain`.
  */
  Vertex_handle insert(const Point & p, Locate_type lt,
                       Face_handle loc, int li, int lj);

  /*!
  Equivalent to `insert(p)`.
  */
  Vertex_handle push_back(const Point& p);

  /*!
  Inserts the points in the iterator range `[first, last)`.
  Returns the number of inserted points. This function uses spatial
  sorting (cf. chapter \ref secspatial_sorting) and therefore is
  not guaranteed to insert the points following the order of
  `InputIterator`. If the third argument `is_large_point_set`
  is set to `true` a heuristic for optimizing the insertion of
  large point sets is applied. \pre The `value_type` of `first` and `last` are `Point`s and all points in the range lie inside `domain`.
  */
  template < class InputIterator > std::ptrdiff_t
  insert(InputIterator first, InputIterator last, bool
         is_large_point_set = false);

  /*!
  Removes the vertex from the triangulation.
  */
  void remove(Vertex_handle v);

/// @}

/// \name Point moving
/// @{

  /*!
  if there is not already another vertex placed on `p`, the
  triangulation is modified such that the new position of vertex
  `v` is `p`, and `v` is returned. Otherwise, the
  triangulation is not modified and the vertex at point `p` is
  returned. \pre `p` lies in the original domain `domain`.
  */
  Vertex_handle move_if_no_collision(Vertex_handle v, const
                                     Point & p);

  /*!
  Moves the point stored in `v` to `p`, while preserving the
  Delaunay property. This performs an action semantically equivalent
  to `remove(v)` followed by `insert(p)`, but is supposedly
  faster to performing these two operations separately when the point
  has not moved much. Returns the handle to the new vertex.
  \pre `p` lies in the original domain `domain`.
  */
  Vertex_handle move_point(Vertex_handle v, const Point &
                           p);

/// @}

/// \name Queries
/// @{

  /*!
  Returns any nearest vertex to the point `p`, or the default constructed
  handle if the triangulation is empty. The optional argument `f` is a hint
  specifying where to start the search. It always returns a vertex
  corresponding to a point inside `domain` even if computing in a
  multiply sheeted covering space.
  \pre `f` is a face of `dt` and `p` lies in the original domain `domain`.

  */
  Vertex_handle nearest_vertex(Point p,
                               Face_handle f = Face_handle());

  /*!
  Returns on which side of the circumcircle of face `f` lies
  the point `p`. The circle is assumed to be counterclockwise
  oriented, so its positive
  side correspond to its bounded side.
  This predicate is available only if the corresponding predicates on
  points is provided in the geometric traits class.
  */
  Oriented_side
  side_of_oriented_circle(Face_handle f, const Point & p);

/// @}

/// \name
/// A point-offset pair (`p`,`off`) is said to be in conflict with a
/// face `f` iff `dt`.`side_of_circle(f, p, off)` returns
/// `ON_BOUNDED_SIDE`. The set of faces that are in conflict with
/// (`p`,`off`) is star-shaped.
/// @{

  /*!
  `OutputItFaces` is an output iterator with `Face_handle` as
  value type. `OutputItBoundaryEdges` stands for an output
  iterator with `Edge` as value type. This members function
  outputs in the container pointed to by `fit` the faces which are
  in conflict with point `p` i. e. the faces whose circumcircle
  contains `p`. It outputs in the container pointed to by
  `eit` the boundary of the zone in conflict with `p`.
  The boundary edges of the conflict zone are output in
  counter-clockwise order and each edge is described through its
  incident face which is not in conflict with `p`. The function
  returns in a std::pair the resulting output
  iterators. \pre `start` is in conflict with `p` and `p` lies in the original domain `domain`.
  */
  template <class OutputItFaces, class OutputItBoundaryEdges>
  std::pair<OutputItFaces, OutputItBoundaryEdges>
  get_conflicts_and_boundary(const Point &p, OutputItFaces fit,
                             OutputItBoundaryEdges eit, Face_handle start) const;

  /*!
  same as above except that only the faces in conflict with `p` are
  output. The function returns the resulting output
  iterator. \pre `start` is in conflict with `p` and `p` lies in the original domain `domain`.
  */
  template <class OutputItFaces> OutputItFaces get_conflicts
  (const Point &p, OutputItFaces fit, Face_handle start) const;

  /*!
  `OutputItBoundaryEdges` stands for an output iterator with
  `Edge` as value type. This function outputs in the container
  pointed to by `eit`, the boundary of the zone in conflict with
  `p`. The boundary edges of the conflict zone are output in
  counterclockwise order and each edge is described through the
  incident face which is not in conflict with `p`. The function
  returns the resulting output iterator. \pre `start` is in conflict with `p` and `p` lies in the original domain `domain`.
  */
  template <class OutputItBoundaryEdges> OutputItBoundaryEdges
  get_boundary_of_conflicts(const Point &p, OutputItBoundaryEdges eit,
                            Face_handle start) const;

/// @}

/// \name Voronoi diagram
/// The following member functions provide the elements of the dual Voronoi diagram.
/// @{

  /*!
  Compute the circumcenter of the face pointed to by f. This function
  is available only if the corresponding function is provided in the
  geometric traits.
  */
  Point circumcenter(Face_handle f) const;

  /*!
  Returns the center of the circle circumscribed to face `f`.
  */
  Point dual(const Face_handle &f) const;

  /*!
  returns a segment whose endpoints are the duals of both incident
  faces.
  */
  Segment dual(const Edge &e) const;

  /*!
  Idem
  */
  Segment dual(const Edge_circulator& ec) const;

  /*!
  Idem
  */
  Segment dual(const Edge_iterator& ei) const;

  /*!
  output the dual Voronoi diagram to stream `ps`.
  */
  template < class Stream> Stream& draw_dual(Stream & ps);

/// @}

/// \name Predicates
/// @{

  /*!
  Returns the side of `p` with respect to the circle circumscribing the
  triangle associated with `f`. Periodic copies are checked if
  necessary.
  */
  Oriented_side side_of_oriented_circle(Face_handle f, const
                                        Point& p, bool perturb ) const;

/// @}

/// \name Checking
/// \cgalAdvancedBegin
/// These methods are mainly a debugging help for the users of advanced features.
/// \cgalAdvancedEnd
/// @{

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  Checks the combinatorial validity of the triangulation and the
  validity of its geometric embedding (see
  Section \ref P2Triangulation2secintro). Also checks that all the
  circumscribing circles of faces are empty.

  When `verbose` is set to true, messages describing the first
  invalidity encountered are printed.
  \cgalAdvancedEnd
  */
  bool
  is_valid(bool verbose = false) const;

  /*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  Checks the combinatorial and geometric validity of the face (see
  Section \ref P2Triangulation2secintro). Also checks that the
  circumscribing circle of faces is empty.

  When `verbose` is set to true, messages are printed to give
  a precise indication of the kind of invalidity encountered.
  \cgalAdvancedEnd
  */
  bool
  is_valid(Face_handle f, bool verbose = false) const;

/// @}

}; /* end Periodic_2_Delaunay_triangulation_2 */
} /* end namespace CGAL */
