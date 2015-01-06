/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// Geometric traits concept for the functions `self_intersect`
class SelfIntersectionTraits{
public:
  /// @name Geometric Types
  /// @{
  /// 3D point type
  typedef unspecified_type Point_3;
  /// 3D triangle type
  typedef unspecified_type Triangle_3;
  /// 3D segment type
  typedef unspecified_type Segment_3;
  /// @}

  /// @name Functors
  /// @{
  /// Functor constructing triangles. It provides `Triangle_3 operator() const(const Point_3&, const Point_3&, const Point_3&)
  typedef unspecified_type Construct_triangle_3;
  /// Functor constructing segments. It provides `Segment_3 operator() const(const Point_3&, const Point_3&)
  typedef unspecified_type Construct_segment_3;
  /// Functor testing intersections between triangles and segment. It provides `bool operator() const (const Triangle_3&, const Segment_3&)` and `bool operator() const (const Triangle_3&, const Triangle_3&)`
  typedef unspecified_type Do_intersect_3;
  /// @}

  /// @name Functions
  /// @{
  Construct_triangle_3 construct_triangle_3_object() const;
  Construct_segment_3 construct_segment_3_object() const;
  Do_intersect_3 do_intersect_3_object() const;
  /// @}
};
