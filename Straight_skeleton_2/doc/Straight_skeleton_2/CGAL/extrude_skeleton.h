namespace CGAL {

/**
* \ingroup PkgStraightSkeleton2Extrusion
*
* \brief constructs the straight skeleton-based extrusion of a polygon with holes.
*
* Given a polygon with holes and a set of weights, the skeleton extrusion is a volume constructed
* from the weighted straight skeleton by associating a height to the vertices of the skeleton,
* which corresponds to the time at the vertex. The input polygon is placed at `z = 0`.
*
* This function allows cropping the extruded skeleton at a maximum height, using the optional
* `maximum_height()` named parameter.
*
* The result is a closed, 2-manifold surface triangle mesh. Note that this mesh can have non-local
* self-intersections if a maximal height is provided due to possible (geometric) non-manifold occurrences.
*
* @tparam PolygonWithHoles must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`)
                           or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).
* @tparam PolygonMesh a model of `MutableFaceGraph`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param polygon the polygon with holes
* @param out the output polygon mesh
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{weights}
*     \cgalParamDescription{a range of weights, one for each edge of the input polygon}
*     \cgalParamType{a model of `Range` whose value type is a model of `Range` whose value type is `FT`}
*     \cgalParamDefault{an empty range (uniform weights are used)}
*     \cgalParamExtra{Weights should be finite and all of the same sign.
*                     Negative weights are used to signify exterior extrusion.
*                     Contrary to weighted skeleton functions, a weight `0` can here be used to signify a vertical extrusion of this edge.}
*     \cgalParamExtra{If neither `weights` nor `angles` are provided, uniform weights are used.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{angles}
*     \cgalParamDescription{a range of angles, one for each edge of the input polygon. For each contour
*                           edge, the angle is the dihedral angle between the `z=0` plane and the desired
*                           extrusion plane.}
*     \cgalParamType{a model of `Range` whose value type is a model of `Range` whose value type is `FT`}
*     \cgalParamDefault{an empty range (uniform weights are used)}
*     \cgalParamExtra{Angles are measured in degrees and should be strictly within `0` and `180` degrees
*                     and should be eitger all acute (inward extrusion) or all obtuse (outward extrusion).}
*     \cgalParamExtra{This parameter is ignored if the `weights` parameter is provided.}
*     \cgalParamExtra{The conversion to weights involves trigonometry and will be inexact,
*                     even when using a number type with exact square roots.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{maximum_height}
*     \cgalParamDescription{the maximum height of the extrusion}
*     \cgalParamType{`FT`, a model of `FieldNumberType` convertible to the kernel type of the polygon.}
*     \cgalParamDefault{unused}
*     \cgalParamExtra{This parameter should not be null, but can be negative; in this case, the polygon is extruded downwards.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{verbose}
*     \cgalParamDescription{Whether information about the extrusion process should be printed to `std::cout`}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`false`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the polygon's point type, using `CGAL::Kernel_traits`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `out`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `InK::Point_2` as value type.}
*     \cgalParamDefault{the internal property map for `CGAL::vertex_point_t` of `out`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* @return `true` if the construction was successful, `false` otherwise.
*
* @pre `polygon` is a weakly simple polygon with holes.
* @pre Holes are weakly simple polygons that do not intersect each other or the outer boundary.
*
* @sa `CGAL::create_interior_weighted_skeleton_and_offset_polygons_with_holes_2()`
* @sa `CGAL::Straight_skeleton_2`
*/
template <typename Polygon,
          typename FT,
          typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
bool extrude_skeleton(const Polygon& polygon,
                      PolygonMesh& out,
                      const NamedParameters& np = parameters::default_values());

} /* namespace CGAL */
