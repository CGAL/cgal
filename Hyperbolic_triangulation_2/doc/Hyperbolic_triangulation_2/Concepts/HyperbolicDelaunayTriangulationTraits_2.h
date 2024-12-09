// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgHyperbolicTriangulation2Concepts
\cgalConcept

\cgalRefines{DelaunayTriangulationTraits_2}

The concept `HyperbolicDelaunayTriangulationTraits_2` describes the set of requirements
to be fulfilled by any class used to instantiate the first template parameter of the class
`CGAL::Hyperbolic_Delaunay_triangulation_2<Traits, Tds>`. It defines the geometric objects
(points, segments...) forming the triangulation together with geometric predicates and
constructions on these objects.

This concept refines `DelaunayTriangulationTraits_2` because the class `CGAL::Hyperbolic_Delaunay_triangulation_2`
internally relies on the class `CGAL::Delaunay_triangulation_2`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Hyperbolic_Delaunay_triangulation_traits_2}
\cgalHasModels{CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2}
\cgalHasModelsEnd
*/


class HyperbolicDelaunayTriangulationTraits_2 {


public:

  /// \name Types
  /// @{

        /*!
                Field number type.
        */
        typedef unspecified_type                          FT;

        /*!
                Represents a point in the Poincaré disk or on the (unit) circle at infinity.
        */
        typedef unspecified_type                         Hyperbolic_point_2;

        /*!
                Represents the dual object of a triangle in the hyperbolic Delaunay triangulation.
                The dual of a Delaunay triangle is the <i>hyperbolic</i> center of the circle circumscribing it.
        */
        typedef unspecified_type                     Hyperbolic_Voronoi_point_2;

        /*!
                Represents a hyperbolic segment defined by two points.
                In the Poincaré disk model, a hyperbolic segment is supported either by the Euclidean circle
                that passes through the two points and is perpendicular to the circle at infinity, or by the
                Euclidean line that passes through the two points and the origin. Abusively, we allow one or
                both endpoints of the segment to lie on the circle at infinity, so a hyperbolic segment can
                actually represent a hyperbolic ray or a hyperbolic line.
        */
        typedef unspecified_type                Hyperbolic_segment_2;

        /*!
                Represents a triangle in the hyperbolic plane defined by three hyperbolic points.
        */
        typedef unspecified_type                         Hyperbolic_triangle_2;
  /// @}

  /// \name Predicate Types
  /// @{

        /*!
                A predicate object. Must provide the function operator

                `Oriented_side operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q, Hyperbolic_point_2 r, Hyperbolic_point_2 t),`
                which returns the position of the point `t` relative to the oriented circle
                defined by the points `p, q`, and `r`.
        */
        typedef unspecified_type                        Side_of_oriented_circle_2;


        /*!
                A predicate object. Must provide the function operator

                `Oriented_side operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q, Hyperbolic_point_2 query),`
                which returns the position of the point `query` relative to the oriented hyperbolic
                segment with vertices `p` and `q`.
        */
        typedef unspecified_type                         Side_of_oriented_hyperbolic_segment_2;


        /*!
                A predicate object. Must provide the function operator

                `bool operator()(Hyperbolic_point_2 p0, Hyperbolic_point_2 p1, Hyperbolic_point_2 p2),`

                which returns a boolean indicating whether the triangle defined
                by the points `p0, p1,` and `p2` is hyperbolic (i.e., if its
                circumscribing disk is contained in the unit disk). It must also
                provide the function operator

                `bool operator() (Hyperbolic_point_2 p0, Hyperbolic_point_2 p1, Hyperbolic_point_2 p2, int& ind),`

                which returns whether the triangle is hyperbolic, and if not stores in `ind` the index of the
                non-hyperbolic edge of the triangle, as defined in \cgalCite{cgal:bdt-hdcvd-14}.
                The edge of the triangle opposite to `pj` for `j = 0,1,2` is considered to have index `j`.

        */
        typedef unspecified_type                         Is_Delaunay_hyperbolic;

  /// @}

  /// \name Construction Types
  /// @{


        /*!
                A constructor object.

                Must provide the function operator

                `Hyperbolic_segment_2 operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q),`

                which constructs a hyperbolic segment from two points `p` and `q`.
                Note that `p` and `q` may also lie on the circle at infinity.
        */
        typedef unspecified_type                         Construct_hyperbolic_segment_2;

        /*!
                A constructor object. Must provide the function operator

                `Hyperbolic_Voronoi_point_2 operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q, Hyperbolic_point_2 r),`

                which constructs the hyperbolic circumcenter of the triangle with
                vertices `p, q`, and `r`.
        */
        typedef unspecified_type                         Construct_hyperbolic_circumcenter_2;

        /*!
                A constructor object. Must provide the function operator

                `Hyperbolic_segment_2 operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q),`

                which constructs the hyperbolic bisector of two points `p` and `q` lying
                in the Poincaré disk. The endpoints of the resulting hyperbolic segment
                lie on the circle at infinity. It must also provide the function operator

                `Hyperbolic_segment_2 operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q, Hyperbolic_point_2 r),`

                where the points `p, q`, and `r` lie in the Poincaré disk. This overloaded
                version constructs the hyperbolic bisector of the segment [p,q] limited by
                the hyperbolic circumcenter of `p, q, r` on one side and the circle at
                infinity on the other. Moreover, it must provide the function operator

                `Hyperbolic_segment_2 operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q, Hyperbolic_point_2 r, Hyperbolic_point_2 s),`

                where the points `p, q, r`, and `s` lie in the Poincaré disk. This overloaded
                version constructs the hyperbolic bisector of the segment [p,q] limited by
             the hyperbolic circumcenter of `p, q, r` on one side, and the hyperbolic
             circumcenter of `p, s, q` on the other side.
        */
        typedef unspecified_type                         Construct_hyperbolic_bisector_2;
  /// @}


  /// \name Operations
  /// The following functions give access to the predicate objects.
  /// @{
        Orientation_2                          orientation_2_object();
        Side_of_oriented_circle_2              side_of_oriented_circle_2_object();
        Side_of_oriented_hyperbolic_segment_2          side_of_oriented_hyperbolic_segment_2_object();
        Is_Delaunay_hyperbolic                                               is_Delaunay_hyperbolic_object();
  /// @}

  /// \name
  /// The following functions must be provided only if the methods of `Hyperbolic_Delaunay_triangulation_2`
  /// that return elements of the Voronoi diagram are instantiated:
  /// @{
        Construct_hyperbolic_segment_2       construct_hyperbolic_segment_2_object();
        Construct_hyperbolic_circumcenter_2  construct_hyperbolic_circumcenter_2_object();
        Construct_hyperbolic_bisector_2      construct_hyperbolic_bisector_2_object();
  /// @}


};

