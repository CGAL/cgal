
/*!
\ingroup PkgFloatSnapRounding2Concepts
\cgalConcept

The concept `FloatSnapRoundingTraits_2` lists the set of requirements that must be fulfilled by
an instance of the `Traits` template-parameter of
the free functions \ref `CGAL::snap_polygons_2()`, `CGAL::compute_snapped_subcurves_2()` and `CGAL::double_snap_rounding_2()`.
The list includes the nested types of the geometric primitives used in this class and
some function object types for the required predicates on those primitives.

\cgalRefines{AosTraits_2}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Float_snap_rounding_traits_2<Kernel>}
\cgalHasModelsEnd
*/

class FloatSnapRoundingTraits_2 {
public:

/// \name Types
/// @{

/*!
The number type. This type must fulfill the requirements on
`FieldNumberType`
*/
typedef unspecified_type FT;

/*!
Models the concept `ArrTraits::Point_2`.
*/
typedef unspecified_type Point_2;

/*!
Models the concept `ArrTraits::XMonotoneCurve_2`.
*/
typedef unspecified_type Segment_2;

/// @}

/// \name Functor Types
/// @{

/*!
Models the concept `Kernel::ConstructSource_2`.
*/
typedef unspecified_type Construct_source_2;

/*!
Models the concept `Kernel::ConstructTarget_2`.
*/
typedef unspecified_type Construct_target_2;

/*!
Models the concept `Kernel::ConstructSegment_2`.
*/
typedef unspecified_type Construct_segment_2;

/*!
Models the concept `Kernel::LessXY_2`.
*/
typedef unspecified_type Less_xy_2;

/*!
Models the concept `Kernel::LessY_2`.
*/
typedef unspecified_type Less_y_2;

/*!
Models the concept `Kernel::Equal_2`.
*/
typedef unspecified_type Equal_2;

/*!
Models the concept `FSRTraits_2::ConstructRoundPoint_2`.
*/
typedef unspecified_type Construct_round_point_2;

/*!
Models the concept `FSRTraits_2::SquaredRoundBound_2`.
*/
typedef unspecified_type Squared_round_bound_2;

/*!
Models the concept `FSRTraits_2::ConverterToExact`.
*/
typedef unspecified_type Converter_to_exact;

/*!
Models the concept `FSRTraits_2::ConverterFromExact`.
*/
typedef unspecified_type Converter_from_exact;

/// @}

/// \name Accessing Functor Objects
/// @{

/*!

*/
Construct_source_2 construct_source_2_object();

/*!

*/
Construct_target_2 construct_target_2_object();

/*!

*/
Construct_segment_2 construct_segment_2_object();

/*!

*/
Less_xy_2 less_xy_2_object();

/*!

*/
Less_y_2 less_y_2_object();

/*!

*/
Construct_round_point_2 construct_round_point_2_object();

/*!

*/
Squared_round_bound_2 squared_round_bound_2_object();


/// @}

}; /* end FloatSnapRoundingTraits_2 */


namespace FSRTraits_2{

/*!
  \ingroup PkgFloatSnapRounding2Concepts
  \cgalConcept
  \cgalHasModelsBegin
  \cgalHasModelsBare{\link FloatSnapRoundingTraits_2::ConstructRoundPoint_2 `Float_snap_rounding_traits_2::Construct_round_point_2` \endlink}
  \cgalHasModelsEnd
*/
class ConstructRoundPoint_2
{
  public:

  /*!
  Given a point, construct its rounded version
  */
  Point_2 operator()(Point_2 p);
};

/*!
  \ingroup PkgFloatSnapRounding2Concepts
  \cgalConcept
  \cgalHasModelsBegin
  \cgalHasModelsBare{\link FloatSnapRoundingTraits_2::SquaredRoundBound_2 `Float_snap_rounding_traits_2::Squared_round_bound_2` \endlink}
  \cgalHasModelsEnd
*/
class SquaredRoundBound_2
{
  public:

  /*!
  Given a point, compute an upper bound of the squared distance between its exact coordinates and its rounded coordinates
  */
  double operator()(Point_2 p);
};

/*!
  \ingroup PkgFloatSnapRounding2Concepts
  \cgalConcept
  \cgalHasModelsBegin
  \cgalHasModelsBare{\link FloatSnapRoundingTraits_2::ConverterToExact `Float_snap_rounding_traits_2::Convertex_to_exact` \endlink}
  \cgalHasModelsEnd
*/
class ConverterToExact
{
  public:
  typedef unspecified_type InputPoint_2;
  typedef unspecified_type InputSegment_2;
  /*!
  Convert input points (i.e. segments) into points (i.e. segments) of the type of the traits
  */
  Point_2 operator()(InputPoint_2 p);
  Segment_2 operator()(InputSegment_2 p);
};

/*!
  \ingroup PkgFloatSnapRounding2Concepts
  \cgalConcept
  \cgalHasModelsBegin
  \cgalHasModelsBare{\link FloatSnapRoundingTraits_2::ConverterToExact `Float_snap_rounding_traits_2::Convertex_to_exact` \endlink}
  \cgalHasModelsEnd
*/
class ConverterFromExact
{
  public:
  typedef unspecified_type OutputPoint_2;
  typedef unspecified_type OutputSegment_2;
  /*!
  Convert points (i.e. segments) of the type of the traits into points (i.e. segments) of the output
  */
  OutputPoint_2 operator()(Point_2 p);
  OutputSegment_2 operator()(Segment_2 p);
};



} /* end of namespace SRTraits_2 */
