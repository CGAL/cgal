namespace CGAL {

/*!
* \ingroup PkgFeatureGraphDetector
*
* Functor for sharp feature detection in labeled images.
*/
struct Detect_sharp_features_in_labelled_image
{
  /*!
  * \brief Detects and constructs the feature graph from a labelled image
  *
  * The feature graph has polylines that lie on the
  * sharp features of the surface defined by the image's subdomains.
  * Each subdomain inside the bounding box
  * of the input labeled image is defined as the set of voxels
  * with the same value. The outside of the bounding box
  * of the image is considered as a subdomain with voxel value `value_outside`.
  *
  * \tparam Point_3 class model of `Kernel::Point_3`.
  * It defines the feature graph point type.
  *
  * \tparam Image class model of `FeatureImage_3`.
  * It represent a labelled image.
  *
  * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * \param image the input image
  *
  * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
  * \cgalNamedParamsBegin
  *   \cgalParamSectionBegin{sharpness_estimator}
  *     \cgalParamDescription{a functor model of `SharpnessEstimator`.
  *       See \ref PkgFeatureGraphSharpnessEstimator for available functors.}
  *     \cgalParamDefault{`CGAL::Sharpness_estimator::AmbrosioTortorelli_on_image`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{selection_threshold}
  *     \cgalParamDescription{a threshold on the sharpness value.
  *       Elements with a sharpness value above this threshold are considered sharp.
  *       The selected elements determines the topology of the output feature graph.}
  *     \cgalParamDefault{`FT(0.25)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{regularization_distance}
  *     \cgalParamDescription{a threshold on the distance between lines.
  *     In the regularization step, if the maximum distance of a line to another line
  *     is less than this threshold, then the line is collapsed.
  *     It can be a constant or a functional.
  *     If it is a functional, it must implement
  *     `FT operator()(const Point_3&) const`.}
  *     \cgalParamDefault{`FT(4.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{do_optimize}
  *     \cgalParamDescription{a boolean indicating if the optimization step is called.}
  *     \cgalParamDefault{`true`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{optimization_parameters}
  *     \cgalParamDescription{an instance of a class model of `OptimizationParameters`.
  *     It describe the parameters for the optimization step.}
  *     \cgalParamDefault{`CGAL::Optimization_parameters_in_image()`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{point_to_element_output_map}
  *     \cgalParamDescription{an output property map that will be filled if supplied.
  *     It allows to retrieve a surface element from the index of a point in the output feature graph.
  *     It must be a model of `WritablePropertyMap`,
  *     so it must implement `put(output_pmap, point_index, std::make_pair(element_type, element_index))`
  *     where the key type is the point index type and the value type is a pair of `CGAL::Surface_element_type` and the element index type.}
  *     \cgalParamDefault{`CGAL::Optimization_parameters_in_image()`}
  *   \cgalParamSectionEnd
  * \cgalNamedParamsEnd
  *
  * \returns a `CGAL::Feature_graph<Point>`
  * containing the constructed features.
  *
  * \sa `Detect_sharp_features_in_surface::oprator()()`
  */
  template<typename Point_3, typename Image, typename CGAL_NP_TEMPLATE_PARAMETERS>
  CGAL::Feature_graph<Point_3>
  operator()(const Image& image, const CGAL_NP_CLASS& np = parameters::default_values()) const;
};

/*!
* \ingroup PkgFeatureGraphDetector
*
* Functor for sharp feature detection in surfaces.
*/
struct Detect_sharp_features_in_surface
{
  /*!
  * \brief Detects and constructs the feature graph from a surface
  *
  * The feature graph has polylines that lie on the
  * sharp features of the surface.
  *
  * \tparam Point_3 a model of `Kernel::Point_3`.
  * It defines the feature graph point type.
  *
  * \tparam Surface a model of `FaceListGraph` that represents a surface mesh.
  *
  * \param surface the input surface
  *
  * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
  * \cgalNamedParamsBegin
  *   \cgalParamSectionBegin{sharpness_measure}
  *     \cgalParamDescription{a functor that implements
  *       `double operator()(const Element& element) const`,
  *       where `Element` is a model of `SurfaceElement_3`,
  *       that returns a sharpness value for that element.
  *       See \ref PkgFeatureGraphSharpnessEstimator for available functors.}
  *     \cgalParamDefault{`parameters::features(domain)`}
  *   \cgalParamSectionEnd
  * \cgalNamedParamsEnd
  * TODO PARAMS
  *
  * \returns a `CGAL::Feature_graph<Point>`
  * containing the constructed features.
  */
  template<typename Point_3, typename Surface, typename CGAL_NP_TEMPLATE_PARAMETERS>
  CGAL::Feature_graph<Point_3>
  operator()(const Surface& surface, const CGAL_NP_CLASS& np = parameters::default_values()) const;
};


} /* namespace CGAL */