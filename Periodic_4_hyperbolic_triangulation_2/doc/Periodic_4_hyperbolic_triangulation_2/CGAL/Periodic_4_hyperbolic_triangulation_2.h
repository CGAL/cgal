// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

namespace CGAL {


/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2MainClasses

The class `Periodic_4_hyperbolic_triangulation_2` offers base functionalities needed by
`Periodic_4_hyperbolic_Delaunay_triangulation_2`. Note that this class does not support
modification of the triangulation (insertion or removal).

The class expects two template parameters.

\tparam GT         Geometric traits parameter. Must be a model of the concept
                        `Periodic_4HyperbolicTriangulationTraits_2`. This parameter has no
                        default value.
\tparam        TDS %Triangulation data structure parameter. Must be a model of the concept
                `TriangulationDataStructure_2` with some additional functionality in the vertices
                and faces, following the concepts `Periodic_4HyperbolicTriangulationVertexBase_2` and
                `Periodic_4HyperbolicTriangulationFaceBase_2`, respectively. The default value for
                this parameter is
                \code
                CGAL::Triangulation_data_structure_2<
                  CGAL::Periodic_4_hyperbolic_triangulation_vertex_base_2<GT>,
                  CGAL::Periodic_4_hyperbolic_triangulation_face_base_2<GT> >
                \endcode
*/
template <         class GT, class TDS        >
class Periodic_4_hyperbolic_triangulation_2 {

public:

        /// \name Types
        /// @{

                typedef GT                                                                                                                 Geometric_traits;
                typedef TDS                                                                                                         Triangulation_data_structure;
        /// @}

        /// \name
        /// The following types represent geometric objects.
        /// @{
                /*!
                        Represents a point in the hyperbolic plane. Note that for Delaunay triangulations
                        of the Bolza surface, all points must lie inside the original hyperbolic octagon.
                */
                typedef typename Triangulation_data_structure::Vertex::Point                 Point;

                /*!
                        Represents a hyperbolic segment in the hyperbolic plane.
                */
                typedef typename Geometric_traits::Hyperbolic_segment_2                         Hyperbolic_segment;

                /*!
                        Represents a triangle contained in the Poincaré disk.
                */
                typedef typename Geometric_traits::Hyperbolic_triangle_2                         Hyperbolic_triangle;
        /// @}

        /// \name
        /// The following type represents hyperbolic translations specific for the Bolza surface.
        /// @{
                typedef typename Geometric_traits::Hyperbolic_translation                         Hyperbolic_translation;
        /// @}

        /// \name
        /// The following types represent images of points, segments, and triangles under
        /// the action of hyperbolic translations.
        /// @{

                /*!
                        Represents a periodic point, i.e., a pair of a point in the original octagon
                        and a hyperbolic translation.
                */
                typedef std::pair<Point, Hyperbolic_translation>                                         Periodic_point;

                /*!
                        Represents a periodic segment, defined by two periodic points.
                */
                typedef std::array< std::pair<Point, Hyperbolic_translation>, 2 >        Periodic_segment;

                /*!
                        Represents a periodic triangle, defined by three periodic points.
                */
                typedef std::array< std::pair<Point, Hyperbolic_translation>, 3 >         Periodic_triangle;
        /// @}


        /// \name
        /// The following types give access to the elements of the triangulation.
        /// @{
                typedef typename Triangulation_data_structure::Edge                                 Edge;
                typedef typename Triangulation_data_structure::Face                                 Face;
                typedef typename Triangulation_data_structure::Vertex_handle                 Vertex_handle;
                typedef typename Triangulation_data_structure::Face_handle                         Face_handle;
        /// @}

        /// \name
        /// The following iterator and circulator types are defined to give access over the vertices,
        /// edges, and faces of the triangulation.
        /// @{
                typedef typename Triangulation_data_structure::Face_iterator        Face_iterator;
                typedef typename Triangulation_data_structure::Edge_iterator        Edge_iterator;
                typedef typename Triangulation_data_structure::Vertex_iterator      Vertex_iterator;
                typedef typename Triangulation_data_structure::Face_circulator      Face_circulator;
                typedef typename Triangulation_data_structure::Edge_circulator      Edge_circulator;
                typedef typename Triangulation_data_structure::Vertex_circulator    Vertex_circulator;
        /// @}

        /// \name
        /// The following enumeration type indicates where a point is located in the triangulation.
        /// @{
        enum Locate_type {
                VERTEX = 0, /*!< the located point coincides with a vertex of the triangulation */
                EDGE,                /*!< the point is in the relative interior of an edge */
                FACE                 /*!< the point is in the interior of a facet */
        };
        /// @}


        /// \name Creation
        /// @{
                /*!
                        %Default constructor, with an optional parameter for the geometric traits object.
                */
                Periodic_4_hyperbolic_triangulation_2( const Geometric_traits &gt = Geometric_traits() );

                /*!
                        Copy constructor.
                */
                Periodic_4_hyperbolic_triangulation_2(const Periodic_4_hyperbolic_triangulation_2& tr);
        /// @}


        /// \name Assignment
        /// @{
                /*!
                        The triangulation `tr` is duplicated, and modifying the copy after the duplication does not modify the original.
                */
                Periodic_4_hyperbolic_triangulation_2& operator=(Periodic_4_hyperbolic_triangulation_2 tr);

                /*!
                        The triangulation is swapped with `tr`.
                */
                void swap(Periodic_4_hyperbolic_triangulation_2& tr);

                /*!
                        Deletes all faces and vertices of the triangulation.
                */
                void clear();

                /*!
                        Equality operator.
                        \todo implement
                */
                bool operator==(const Periodic_4_hyperbolic_triangulation_2<GT, TDS>& tr1,
                                                const Periodic_4_hyperbolic_triangulation_2<GT, TDS>& tr2);

                /*!
                        Inequality operator.
                        \todo implement
                */
                bool operator!=(const Periodic_4_hyperbolic_triangulation_2<GT, TDS>& tr1,
                                                const Periodic_4_hyperbolic_triangulation_2<GT, TDS>& tr2);
        /// @}



        /// \name Access functions
        /// @{

                /*!
                        Returns a const reference to the geometric traits object.
                */
                const Geometric_traits& geom_traits() const;

                /*!
                        Returns a const reference to the triangulation data structure.
                */
                const Triangulation_data_structure & tds() const;

                /*!
                        Returns a reference to the triangulation data structure.
                */
                Triangulation_data_structure & tds();

                  /*!
                        Returns the number of vertices in the triangulation.
                */
                size_type number_of_vertices() const;

                /*!
                        Returns the number of edges in the triangulation.
                */
                size_type number_of_edges() const;

                /*!
                        Returns the number of faces in the triangulation.
                */
                size_type number_of_faces() const;
        /// @}



        /// \name Geometric access functions
        /// The functions below return periodic points, segments, and triangles.
        /// @{
                /*!
                        Returns the periodic point given by the `i`-th vertex of face `f`, that is the point
                        in the original domain and the translation of the vertex in `f`.
                        \pre \f$0 \leq i \leq 2\f$
                */
                Periodic_point periodic_point( const Face_handle f, int i) const;

                /*!
                        Returns the periodic segment formed by the two point-translation pairs `(p1, tr1)` and `(p2, tr2)`.
                */
                Periodic_segment periodic_segment(        const Point& p1,                                         const Point& p2,
                                                                                        const Hyperbolic_translation& tr1,         const Hyperbolic_translation& tr2) const;

                /*!
                        Returns the periodic segment formed by the two point-translation pairs `(p1, Id)` and `(p2, Id)`,
                        where `Id` is the identity translation.
                */
                Periodic_segment periodic_segment(        const Point& p1,                                         const Point& p2) const;

                /*!
                        Returns the periodic segment formed by the endpoints of edge `(f,i,j)`.
                        \pre \f$ 0 \leq i,j \leq 2, \qquad i \neq j \f$
                */
                Periodic_segment periodic_segment(const Face_handle c, int i, int j) const;

                  /*!
                        Returns the periodic segment formed by the endpoints of edge `e`.
                        Note that the translations in the resulting periodic segment are determined
                        by `e.first`.
                  */
                Periodic_segment periodic_segment(const Edge & e) const;

                /*!
                        Returns the periodic segment formed by the endpoints of edge `e` translated by `tr`.
                        Note that the translations in the resulting segment are determined by the translations
                        in `e.first`, multiplied on the left by `tr`.
                */
                Periodic_segment periodic_segment(const Edge & e, const Hyperbolic_translation& tr) const;

                /*!
                        Returns the periodic triangle formed by the pairs of points and translations
                        corresponding to each vertex of `f`.
                */
                Periodic_triangle periodic_triangle(const Face & f) const;

                /*!
                        Returns the periodic triangle formed by the pairs of points and translations
                        corresponding to each vertex of `f`, translated by `tr`.
                */
                Periodic_triangle periodic_triangle(const Face & f, const Hyperbolic_translation & tr) const;
        /// @}


        /// \name
        /// The functions below return non-periodic points, segments and triangles.
        /// @{
                /*!
                        Converts the periodic point `pp` into a point in the hyperbolic plane by applying `pp.second` to `pp.first`.
                */
                Point construct_point(const Periodic_point & pp) const;

                /*!
                        Constructs the hyperbolic segment formed by the endpoints of edge `(fh, idx)`.
                        \pre \f$ 0 \leq idx \leq 2 \f$
                */
                Hyperbolic_segment construct_hyperbolic_segment(const Face_handle & fh, int idx) const;

                /*!
                        Returns the hyperbolic segment with endpoints `p1` and `p2`.
                */
                Hyperbolic_segment construct_hyperbolic_segment(const Point& p1, const Point& p2) const;

                /*!
                        Returns the hyperbolic segment formed by the endpoints of `e`.
                */
                Hyperbolic_segment construct_hyperbolic_segment(const pair<Face_handle, int> & e) const;

                /*!
                        Returns the hyperbolic segment formed by the endpoints of `ps`.
                */
                Hyperbolic_segment construct_hyperbolic_segment(const Periodic_segment & ps) const;

                /*!
                        Returns the triangle formed by the vertices of `fh`.
                */
                Triangle construct_triangle(const Face_handle & fh) const;

                /*!
                        Returns the triangle formed by the vertices of `pt`.
                */
                Triangle construct_triangle(const Periodic_triangle & pt) const;
        /// @}




        /// \name Point location
        /// @{

                /*!
                        Returns the face `rf` for which the periodic triangle `(rf, lo)` contains the query
                        point `p`. The parameter `start`, if provided, is used as a starting point for the location.
                */
                Face_handle hyperbolic_periodic_locate(        const Point& p,
                                                                                                Hyperbolic_translation& lo,
                                                                                                const Face_handle start = Face_handle()) const;

                /*!
                        Same as above. The value of the variable `lt` indicates whether        `p` has been located
                        inside the hyperbolic triangle, on one of its sides, or on one of its vertices. If `p`
                        is located on a side or on a vertex, then `li` contains the index of the corresponding
                        edge or vertex in the returned face. The parameter `start`, if provided, is used as a
                        starting point for the location.
                */
                Face_handle hyperbolic_periodic_locate(        const Point& p,
                                                                                                Locate_type& lt,
                                                                                                int& li,
                                                                                                Hyperbolic_translation& lo,
                                                                                                const Face_handle start = Face_handle()) const;

                /*!
                        Returns the canonical representative of the face that contains the query point `p`.
                        The parameter `start`, if provided, is used as a starting point for the location.
                */
                Face_handle hyperbolic_locate(const Point& p,
                                                                          Face_handle start = Face_handle()) const;

                /*!
                        Same as above. The value of the variable `lt` indicates whether        `p` has been located
                        inside the canonical representative, on one of its sides, or on one of its vertices.
                        If `p` is located on a side or on a vertex, then `li` contains the index of the
                        corresponding edge or vertex in the returned face. The parameter `start`, if provided,
                        is used as a starting point for the location.
                */
                Face_handle hyperbolic_locate(const Point& p,
                                                                          Locate_type& lt,
                                                                          int& li,
                                                                          Face_handle start = Face_handle()) const;
        /// @}



        /// \name Predicate functions
        /// @{

                /*!
                        Returns the `Orientation` of the points `p1, p2, p3`.
                        \sa CGAL::orientation()
                */
                Orientation orientation(const Point &p1, const Point &p2, const Point &p3) const;

                /*!
                        Returns the `Orientation` of the three periodic points `(p1, tr1), (p2, tr2)` and `(p3, tr3)`.
                */
                Orientation orientation(const Point  &p1, const Point  &p2, const Point  &p3,
                                                                const Hyperbolic_translation &tr1,
                                                                const Hyperbolic_translation &tr2,
                                                                const Hyperbolic_translation &tr3) const;

                /*!
                        Returns the `Oriented_side` on which `q` is located with respect to the circle defined by the points `p1,p2,p3`.
                        \sa CGAL::side_of_oriented_circle()
                */
                Oriented_side side_of_oriented_circle(        const Point  &p1, const Point &p2,
                                                                                                const Point  &p3, const Point &q) const;

                /*!
                        Returns the `Oriented_side` on which the periodic point `(q, trq)` is located
                        with respect to the circle defined by the periodic points `(p1, tr1), (p2, tr2)`
                        and `(p3, tr3)`.
                */
                Oriented_side side_of_oriented_circle(        const Point  &p1,                                         const Point  &p2,
                                                                                                const Point  &p3,                                         const Point &q,
                                                                                                  const Hyperbolic_translation &tr1,         const Hyperbolic_translation &tr2,
                                                                                                  const Hyperbolic_translation &tr3,         const Hyperbolic_translation &trq) const;
        /// @}



        /// \name Queries
        /// @{

                /*!
                        Tests whether `v` is a vertex of the triangulation
                */
                bool is_vertex(Vertex_handle v) const;

                /*!
                        Tests whether `(u, v)` is an edge of the triangulation.
                        If the edge is found, then it is the `i`-th edge of `fh`.
                */
                bool is_edge(Vertex_handle u, Vertex_handle v,
                                         Face_handle & fh, int & i) const;

                /*!
                        Tests whether there exists a face in the triangulation with vertices
                        `u, v` and `w`. If the face is found, then it is returned as `fh`.
                */
                bool is_face(        Vertex_handle u, Vertex_handle v, Vertex_handle w,
                                                Face_handle & fh) const;

                /*!
                        Tests whether face `f` has `v` as a vertex. If the answer is `true`,
                        then `i` contains the index of `v` in `f`.
                */
                bool has_vertex(const Face_handle f, const Vertex_handle v, const int i) const;
        /// @}



        /// \name Vertex, Edge, and Face iterators
        /// @{

                /*!
                        Starts at an arbitrary vertex. Iterates over all the vertices in the triangulation.
                */
                Vertex_iterator vertices_begin() const;

                /*!
                        Past-the-end iterator.
                */
                Vertex_iterator vertices_end() const;

                /*!
                        Starts at an arbitrary edge. Iterates over all the edges in the triangulation.
                */
                Edge_iterator edges_begin() const;

                /*!
                        Past-the-end iterator.
                */
                Edge_iterator edges_end() const;

                /*!
                        Starts at an arbitrary face. Iterates over all the faces of the triangulation.
                */
                Face_iterator faces_begin() const;

                /*!
                        Past-the-end iterator.
                */
                Face_iterator faces_end() const;
        /// @}


        /// \name Vertex, Edge and Face circulators
        /// @{
                /*!
                        Starts at an arbitrary vertex incident to `v`.
                */
                Vertex_circulator adjacent_vertices(Vertex_handle v) const;

                /*!
                        Starts at the first vertex of `f` adjacent to `v` in counterclockwise order around `v`.
                */
                Vertex_circulator adjacent_vertices(Vertex_handle v, Face_handle f) const;

                /*!
                        Starts at an arbitrary edge incident to `v`.
                */
                Edge_circulator incident_edges(Vertex_handle v) const;

                /*!
                        Starts at the first edge of `f` incident to `v`, in counterclockwise order around `v`.
                */
                Edge_circulator incident_edges(Vertex_handle v, Face_handle f) const;

                /*!
                        Starts at an arbitrary face incident to `v`.
                */
                Face_circulator incident_faces(Vertex_handle v) const;

                /*!
                        Starts at face `f`.
                */
                Face_circulator incident_faces(Vertex_handle v, Face_handle f) const;
        /// @}


        /// \name Traversal of the Vertices, Edges and Faces incident to a vertex
        /// @{

                /*!
                        Copies the faces incident to `v` into the output iterator `faces`.
                */
                template <class OutputIterator>
                OutputIterator incident_faces(Vertex_handle v, OutputIterator faces) const;

                /*!
                        Copies the edges incident to `v` into the output iterator `edges`.
                */
                template <class OutputIterator>
                OutputIterator incident_edges(Vertex_handle v, OutputIterator edges) const;

                /*!
                        Copies the vertices incident to `v` into the output iterator `vertices`.
                */
                template <class OutputIterator>
                OutputIterator adjacent_vertices(Vertex_handle v, OutputIterator vertices) const;

                /*!
                        Returns the degree of `v`, i.e., the number of edges incident to `v`.
                */
                size_type degree(Vertex_handle v) const;
        /// @}



        /// \name Traversal between adjacent faces
        /// @{

                /*!
                        Returns the index of `f` in its `i`-th neighbor.
                */
                int mirror_index(Face_handle f, int i) const;

                /*!
                        Returns the vertex of the `i`-th neighbor of `f` that is opposite to `f`.
                */
                Vertex_handle mirror_vertex(Face_handle f, int i) const;

                /*!
                        Returns the same edge seen from the other adjacent face.
                */
                Edge mirror_edge(Edge e) const;

                /*!
                        Returns the hyperbolic translation for which the `i`-th neighbor of `fh`
                        is adjacent to (the canonical representative of) `fh` in the hyperbolic plane.
                */
                Hyperbolic_translation neighbor_translation(const Face_handle fh, int i) const;
        /// @}


        /// \name Validity checking
        /// @{

                /*!
                        Checks the combinatorial validity of the triangulation, and also the
                        validity of its geometric embedding.
                */
                bool is_valid(bool verbose = false) const;

                /*!
                        Checks the combinatorial validity of face `f`.
                */
                bool is_valid(Face_handle f, bool verbose = false) const;
        /// @}

}; // class Periodic_4_hyperbolic_triangulation_2





}  // namespace CGAL






