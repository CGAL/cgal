/// \ingroup PkgPolygonMeshProcessingConcepts
/// \cgalConcept
///
/// The concept `PMPDistanceTraits` is a refinement of the
/// concepts `AABBGeomTraits_3` and `SpatialSortingTraits_3`. In addition to the types required by
/// those concepts,  it also requires types and functors needed by the functions `approximate_max_distance_to_point_set()`,
/// `sample_triangle_mesh()`, `approximate_Hausdorff_distance()` and `max_distance_to_triangle_mesh()`
///
/// \cgalRefines{AABBGeomTraits_3, SpatialSortingTraits_3}
/// \cgalHasModelsBegin
/// \cgalHasModelsBare{Any 3D Kernel is a model of this concept}
/// \cgalHasModelsEnd

class PMPDistanceTraits{
public:
    /*!
     * A number type model of `Field` and `RealEmbeddable`
     */
    typedef unspecified_type FT;
    ///
    /*! 3D point type
     * It must be default constructible, and can be constructed from 3 objects of type `FT`.
     * `bool operator<(Point_3, Point_3)` to lexicographically compare two points must be available.
     * Access to %Cartesian coordinates must be possible using `Point_3::x()`, `Point_3::y(), Point_3::z()` and
     * `FT operator[](int i)` with  `0 <= i < 3`.
     *
     * There must be a specialization of `CGAL::Kernel_traits` such that
     * `CGAL::Kernel_traits<Point_3>::%Kernel` is a model implementing this concept.
     */
    typedef unspecified_type Point_3;

    /// 3D vector type
    typedef unspecified_type Vector_3;

    /// @name Functors
    /// @{

    /// Functor for computing squared area of a triangle.
    /// It provides `FT operator()(const Point_3&, const Point_3&, const Point_3&) const`
    /// and `FT operator()(const Triangle_3&) const` and has `FT` as result_type.
    typedef unspecified_type Compute_squared_area_3;
    /// Functor for computing squared length of a segment.
    /// and `FT operator()(const Segment_3&) const` and has `FT` as result_type.
    typedef unspecified_type Compute_squared_length_3;
    /// Functor for constructing translated points.
    /// It provides `Point_3 operator()(const Point_3 &, const Vector_3 &)`
    typedef unspecified_type Construct_translated_point_3;
    /// Functor for constructing vectors.
    /// It provides `Vector_3 operator()(const Point_3 &, const Point_3 &)`
    typedef unspecified_type Construct_vector_3;
    /// Functor for constructing scaled vectors.
    /// It provides `Vector_3 operator()(const Vector_3 &, const FT &)`
    typedef unspecified_type Construct_scaled_vector_3;
    /// @}

    /// @name Functions
    /// @{
    Compute_squared_area_3 compute_squared_area_3_object();
    Compute_squared_length_3 compute_squared_length_3_object();
    Construct_translated_point_3 construct_translated_point_3_object();
    Construct_vector_3 construct_vector_3_object();
    Construct_scaled_vector_3 construct_scaled_vector_3_object();
    /// @}
};
