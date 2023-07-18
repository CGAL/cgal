/*!
  \ingroup PkgOrthtreeConcepts
  \cgalConcept

  The concept `OrthtreeTraits` defines the requirements for the
  template parameter of the `CGAL::Orthtree` class.

  \cgalHasModelsBegin
  \cgalHasModels{CGAL::Orthtree_traits_2<GeomTraits>}
  \cgalHasModels{CGAL::Orthtree_traits_3<GeomTraits>}
  \cgalHasModels{CGAL::Orthtree_traits_d<GeomTraits,Dimension>}
  \cgalHasModelsEnd
*/
class OrthtreeTraits
{
public:

  /// \name Types
  /// @{

  typedef unspecified_type Dimension; ///< Dimension type (see `CGAL::Dimension_tag`).
  typedef unspecified_type Bbox_d; ///< Bounding box type.
  typedef unspecified_type FT; ///< The number type of the %Cartesian coordinates of types `Point_d`
  typedef unspecified_type Point_d; ///< Point type.
  typedef unspecified_type Sphere_d; ///< The sphere type for neighbor queries.

  /*!
    A random access iterator type to enumerate the
    %Cartesian coordinates of a point.
  */
  typedef unspecified_type Cartesian_const_iterator_d;
  typedef std::array<FT, Dimension::value> Array; ///< Array used for easy point constructions.

  typedef unspecified_type Adjacency; ///< Specify the adjacency directions

  /*!
    Functor with an operator to construct a `Point_d` from an `Array` object.
  */
  typedef unspecified_type Construct_point_d_from_array;

  /*!
    Functor with an operator to construct a `Bbox_d` from two `Array` objects (coordinates of minimum and maximum points).
  */
  typedef unspecified_type Construct_bbox_d;

  /// @}

  /// \name Operations
  /// @{

  /*!
    Function used to construct an object of type `Construct_point_d_from_array`.
  */
  Construct_point_d_from_array construct_point_d_from_array_object() const;

  /*!
    Function used to construct an object of type `Construct_bbox_d`.
  */
  Construct_bbox_d construct_bbox_d_object() const;

  /// @}
};
