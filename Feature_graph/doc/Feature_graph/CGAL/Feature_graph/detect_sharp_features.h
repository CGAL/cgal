namespace CGAL {

namespace Feature_graph {

/*!
* \ingroup PkgFeatureGraphDetector
*
* Functor for sharp feature detection in labeled images.
*
* \sa `Detect_sharp_features_on_surface`
*/
struct Detect_sharp_features_on_labeled_image
{
  /*!
  * \brief detects and constructs the feature graph from a labeled image
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
  * It represent a labeled image.
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
  *     \cgalParamDefault{`CGAL::Feature_graph::AmbrosioTortorelli_on_image(image).sharpness_functor()`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{regularization_parameters}
  *     \cgalParamDescription{an instance of `CGAL::Feature_graph::Regularization_parameters`.
  *       It describes the parameters for the regularization step.}
  *     \cgalParamDefault{`CGAL::Feature_graph::Regularization_parameters_on_image()`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{do_optimize}
  *     \cgalParamDescription{a boolean indicating if the optimization step is called.}
  *     \cgalParamDefault{`true`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{optimization_parameters}
  *     \cgalParamDescription{an instance of `CGAL::Feature_graph::Optimization_parameters`.
  *       It describe the parameters for the optimization step.}
  *     \cgalParamDefault{`CGAL::Feature_graph::Optimization_parameters_on_image<>()`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{point_to_element_output_map}
  *     \cgalParamDescription{an output property map that will be filled if supplied.
    *     It allows retrieving of the surface element where a feature graph point is embedded.
    *     It must be a model of `WritablePropertyMap`,
    *     so it must implement `put(output_pmap, point_index, element_index)`
    *     where the key type is the point index type of the feature graph, and the value type the element index type of the surface.
    *     The element index correspond to a facet element with the tag `CGAL::DimensionTag<2>`}
  *     \cgalParamDefault{`parameters::default()`}
  *   \cgalParamSectionEnd
  * \cgalNamedParamsEnd
  *
  * \returns a graph model of `VertexAndEdgeListGraph`
  * containing the constructed features.
  */
  template<typename Point_3, typename Image, typename CGAL_NP_TEMPLATE_PARAMETERS>
  unspecified_type operator()(const Image& image, const CGAL_NP_CLASS& np = parameters::default_values()) const;
};

/*!
* \ingroup PkgFeatureGraphDetector
*
* Functor for sharp feature detection in surfaces.
*
* \sa `Detect_sharp_features_on_labeled_image`
*/
struct Detect_sharp_features_on_surface
{
  /*!
  * \brief detects and constructs the feature graph from a surface
  *
  * The feature graph has polylines that lie on the
  * sharp features of the surface.
  *
  * \tparam Point_3 a model of `Kernel::Point_3`.
  * It defines the feature graph point type.
  *
  * \tparam PolygonMesh a model of `FaceListGraph` that represents a surface mesh.
  *
  * \param pmesh a polygon mesh from where the sharp features are extracted.
  *
  * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
  * \cgalNamedParamsBegin
  *   \cgalParamSectionBegin{sharpness_estimator}
  *     \cgalParamDescription{a functor model of `SharpnessEstimator`.
  *       See \ref PkgFeatureGraphSharpnessEstimator for available functors.}
  *     \cgalParamDefault{`CGAL::Feature_graph::Sharpness_estimator::Sharpness_estimator_on_surface(surface)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{regularization_parameters}
  *     \cgalParamDescription{an instance of `CGAL::Feature_graph::Regularization_parameters`.
  *       It describes the parameters for the regularization step.}
  *     \cgalParamDefault{`CGAL::Feature_graph::Regularization_parameters_on_surface()`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{do_optimize}
  *     \cgalParamDescription{a boolean indicating if the optimization step is called.}
  *     \cgalParamDefault{`true`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{optimization_parameters}
  *     \cgalParamDescription{an instance of `Optimization_parameters`.
  *     It describe the parameters for the optimization step.}
  *     \cgalParamDefault{`CGAL::Feature_graph::Optimization_parameters_on_surface<>()`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{point_to_element_output_map}
  *     \cgalParamDescription{an output property map that will be filled if supplied.
  *     It allows retrieving of the surface element where a feature graph point is embedded.
  *     It must be a model of `WritablePropertyMap`,
  *     so it must implement `put(output_pmap, point_index, element_index)`
    *     where the key type is the point index type of the feature graph, and the value type the element index type of the surface.
    *     The element index correspond to a facet element with the tag `CGAL::DimensionTag<2>`}
  *     \cgalParamDefault{`parameters::default()`}
  *   \cgalParamSectionEnd
  * \cgalNamedParamsEnd
  *
  * \returns a graph model of `VertexAndEdgeListGraph`
  * containing the constructed features.
  */
  template<typename Point_3, typename PolygonMesh, typename CGAL_NP_TEMPLATE_PARAMETERS>
  unspecified_type operator()(const PolygonMesh& pmesh, const CGAL_NP_CLASS& np = parameters::default_values()) const;
};

} /* namespace Feature_graph */

} /* namespace CGAL */