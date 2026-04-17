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
  *     \cgalParamDefault{`CGAL::AmbrosioTortorelli_on_image(image).get_sharpness_functor()`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{regularization_distance}
  *     \cgalParamDescription{a threshold on the distance between lines.
  *     In the regularization step, if the maximum distance of a line to another line
  *     is less than this threshold, then the line is collapsed.
  *     It can be a constant or a functor.
  *     If it is a functor, it must implement
  *     <UL>
  *       <LI> `template <typename Element_type_tag, typename Index>`
  *       <BR> `FT operator()(const Point_3& point_in_space, const Index& element_index) const`
  *     </UL>}
  *     \cgalParamDefault{`FT(4.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{do_optimize}
  *     \cgalParamDescription{a boolean indicating if the optimization step is called.}
  *     \cgalParamDefault{`true`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{optimization_parameters}
  *     \cgalParamDescription{an instance of `Optimization_parameters`.
  *     It describe the parameters for the optimization step.}
  *     \cgalParamDefault{`CGAL::Optimization_parameters_in_image()`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{point_to_element_output_map}
  *     \cgalParamDescription{an output property map that will be filled if supplied.
  *     It allows to retrieve the surface element where a feature graph point is embedded.
  *     It must be a model of `WritablePropertyMap`,
  *     so it must implement `put(output_pmap, point_index, element_index)`
  *     where the key type is the point index type and the value type the element index type.
  *     The element index correspond to surface element with the tag `CGAL::Element_type::Surface`}
  *     \cgalParamDefault{`parameters::default()`}
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
  *   \cgalParamSectionBegin{sharpness_estimator}
  *     \cgalParamDescription{a functor model of `SharpnessEstimator`.
  *       See \ref PkgFeatureGraphSharpnessEstimator for available functors.}
  *     \cgalParamDefault{`CGAL::Sharpness_estimator::Sharpness_estimator_on_surface(surface)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{regularization_distance}
  *     \cgalParamDescription{a threshold on the distance between lines.
  *     In the regularization step, if the maximum distance of a line to another line
  *     is less than this threshold, then the line is collapsed.
  *     It can be a constant or a functional.
  *     If it is a functor, it must implement
  *     <UL>
  *       <LI> `template <typename Element_type_tag, typename Index>`
  *       <BR> `FT operator()(const Point_3& point_in_space, const Index& element_index) const`
  *     </UL>}
  *     \cgalParamDefault{`FT(0.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{do_optimize}
  *     \cgalParamDescription{a boolean indicating if the optimization step is called.}
  *     \cgalParamDefault{`true`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{optimization_parameters}
  *     \cgalParamDescription{an instance of `Optimization_parameters`.
  *     It describe the parameters for the optimization step.}
  *     \cgalParamDefault{`CGAL::Optimization_parameters_in_surface()`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{point_to_element_output_map}
  *     \cgalParamDescription{an output property map that will be filled if supplied.
  *     It allows to retrieve the surface element where a feature graph point is embedded.
  *     It must be a model of `WritablePropertyMap`,
  *     so it must implement `put(output_pmap, point_index, element_index)`
  *     where the key type is the point index type and the value type the element index type.
  *     The element index correspond to surface element with the tag `CGAL::Element_type::Surface`}
  *     \cgalParamDefault{`parameters::default()`}
  *   \cgalParamSectionEnd
  * \cgalNamedParamsEnd
  *
  * \returns a `CGAL::Feature_graph<Point>`
  * containing the constructed features.
  */
  template<typename Point_3, typename Surface, typename CGAL_NP_TEMPLATE_PARAMETERS>
  CGAL::Feature_graph<Point_3>
  operator()(const Surface& surface, const CGAL_NP_CLASS& np = parameters::default_values()) const;
};


} /* namespace CGAL */