namespace CGAL {

namespace Feature_graph {

/*!
* \ingroup PkgFeatureGraphNormalEstimator, PkgFeatureGraphSharpnessEstimator
*
* \brief Class that evaluates the Ambriosio-Tortorelli normals and sharpness measure from an image.
* Two estimators can then be retrieved to access the normal and sharpness estimations.
*
* \tparam Vector_3 the type of the normal vector model of `Kernel::Vector_3`.
*
*/
template <typename Vector_3>
struct Image_AmbrosioTortorelli
{
public:
  /// \name Types
  /// @{

  /*!
  * The type of the normal vector.
  */
  typedef Vector_3 Normal_type;

  /*!
  * The type of the estimator that allow the user to retrieve the sharpness values.
  * \cgalModels{SharpnessEstimator}
  */
  typedef unspecified_type Sharpness_estimator;
  /*!
  * The type of the estimator that allow the user to retrieve the normals.
  * \cgalModels{NormalEstimator}
  */
  typedef unspecified_type Normal_estimator;

  /// @}

  /// \name Constructor
  /// @{

  /*!
  * evaluates the normal and sharpness values using the Ambrosio-Tortorelli energy optimization.
  *
  * \param image the image domain.
  * \param selection_threshold a threshold on the sharpness value.
  *     Elements with a sharpness value lower than this threshold are considered flat
  *     and will be given a negative value.
  */
  Image_AmbrosioTortorelli(const CGAL::Image_3& image, const Sharpness_estimator::Sharpness_number_type& selection_threshold = Sharpness_estimator::Sharpness_number_type(0.25));

  /// @}

  /// \name Estimators
  /// @{

  /*!
  * returns the estimator that allow the user to retrieve the sharpness values.
  */
  Sharpness_estimator sharpness_estimator() const;
  /*!
  * returns the estimator that allow the user to retrieve the normals.
  */
  Normal_estimator normal_estimator() const;

  /// @}
};

} /* namespace Feature_graph */

} /* namespace CGAL */