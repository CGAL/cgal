namespace CGAL {

namespace Feature_graph {

/*!
* \ingroup PkgFeatureGraphDetector
*
* Detector for sharp feature extraction in labeled images.
* It constructs a feature graph model of `VertexAndEdgeListGraph`.
*
* \tparam K the requirements for the geometric objects, model of `Kernel`.
*
* \sa `Detect_sharp_features_on_surface`
*/
template <typename K>
struct Detect_sharp_features_on_labeled_image
{
  /// \name Types
  /// @{

  /*!
  * Output feature graph point type
  */
  typedef typename K::Point_3 Point_3;

  /*!
  * Natural number type.
  */
  typedef std::size_t Size;

  /*!
  * Numerical type.
  */
  typedef typename K::FT FT;

  /// @}

  /// \name Detector
  /// @{

  /*!
  * \brief detects and constructs the feature graph from a labeled image
  *
  * The feature graph has polylines that lie on the
  * sharp features of the surface defined by the image's subdomains.
  * Each subdomain inside the bounding box
  * of the input labeled image is defined as the set of voxels
  * with the same value. The outside of the bounding box
  * of the image is considered as a subdomain with voxel value
  * specified by `parameters::value_outside`.
  *
  * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * \param image the input image domain that represent a labeled image.
  *
  * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
  * \cgalNamedParamsBegin
  *   \cgalParamSectionBegin{value_outside}
  *      \cgalParamDescription{the value attached to voxels
  *                            outside the domain where feature are extracted.}
  *      \cgalParamDefault{0}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{sharpness_estimator}
  *     \cgalParamDescription{an estimator model of `SharpnessEstimator`.
  *       See \ref PkgFeatureGraphSharpnessEstimator for available estimators.}
  *     \cgalParamDefault{`CGAL::Feature_graph::Image_AmbrosioTortorelli(image).sharpness_estimator()`}
  *   \cgalParamSectionEnd
  *
  *   \cgalParamSectionBegin{line_distance_threshold}
  *     \cgalParamDescription{the threshold on the distance between lines.
  *         In the regularization step, if the maximum distance of a line to another line
  *         is less than this threshold, then the line is collapsed.
  *         It can be a constant or a scalar field that return the threshold at a point on the surface.
  *         If it is a scalar field, it must implement
  *         <UL>
  *           <LI> `template <typename DimensionTag>`
  *           <BR> `FT operator()(const Point_3& point_in_space, const std::size_t& element_descriptor, const CGAL::Image_3& image) const`
  *         </UL>
  *         <LI> The type `DimensionTag` represents the element type, it can either be `CGAL::Dimension_tag<0>`, `CGAL::Dimension_tag<1>` or `CGAL::Dimension_tag<2>`,
  *         which represent either a vertex element, an edge element, or a face element respectively.
  *         <LI> The parameter `element_descriptor` indicates the element of the image :
  *         from the centroid `(x,y,z)` of the point, line or surface element in the image space,
  *         and given the dimension of the image `(xdim, ydim, zdim)`, the descriptor is
  *         `2(z*xdim*ydim+y*xdim+x)`}
  *     \cgalParamDefault{`FT(4.0)`}
  *     \cgalParamExtra{This parameter sets a baseline distance and is overridden
  *                     by the parameters `parameters::isthmus_line_distance_threshold`, `parameters::simple_line_distance_threshold` and `parameters::corner_line_distance_threshold`
  *                     if they are set.}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{isthmus_line_distance_threshold}
  *     \cgalParamDescription{the threshold on the distance between an isthmus line and any line.
  *         An isthmus line is bounded by a corner that is incident to only the same line.
  *         In the regularization step, if the maximum distance of a line to another line
  *         is less than this threshold, then the line is collapsed.
  *         It can be a constant or a scalar field that return the threshold at a point on the surface.
  *         If it is a scalar field, it must implement
  *         <UL>
  *           <LI> `template <typename DimensionTag>`
  *           <BR> `FT operator()(const Point_3& point_in_space, const std::size_t& element_descriptor, const CGAL::Image_3& image) const`
  *         </UL>
  *         See the description of `parameters::line_distance_threshold` for more details on this scalar field.}
  *     \cgalParamDefault{`FT(-1.0)`}
  *     \cgalParamExtra{Negative values implies the use of the parameter `parameters::line_distance_threshold`.}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{simple_line_distance_threshold}
  *     \cgalParamDescription{the threshold on the distance between a simple line and any line.
  *         A simple line is a line that can be removed without splitting the graph of its neighbors.
  *         In the regularization step, if the maximum distance of a line to another line
  *         is less than this threshold, then the line is collapsed.
  *         It can be a constant or a scalar field that return the threshold at a point on the surface.
  *         If it is a scalar field, it must implement
  *         <UL>
  *           <LI> `template <typename DimensionTag>`
  *           <BR> `FT operator()(const Point_3& point_in_space, const std::size_t& element_descriptor, const CGAL::Image_3& image) const`
  *         </UL>
  *         See the description of `parameters::line_distance_threshold` for more details on this scalar field.}
  *     \cgalParamDefault{`FT(-1.0)`}
  *     \cgalParamExtra{Negative values implies the use of the parameter `parameters::line_distance_threshold`.}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{corner_line_distance_threshold}
  *     \cgalParamDescription{the threshold on the distance between a corner line and any line.
  *         A corner line is neither isthmus nor simple.
  *         In the regularization step, if the maximum distance of a line to another line
  *         is less than this threshold, then the line is collapsed.
  *         It can be a constant or a scalar field that return the threshold at a point on the surface.
  *         If it is a scalar field, it must implement
  *         <UL>
  *           <LI> `template <typename DimensionTag, typename Descriptor>`
  *           <BR> `FT operator()(const Point_3& point_in_space, const Descriptor& element_descriptor, const CGAL::Image_3& image) const`
  *         </UL>
  *         See the description of `parameters::line_distance_threshold` for more details on this scalar field.}
  *     \cgalParamDefault{`FT(-1.0)`}
  *     \cgalParamExtra{Negative values implies the use of the parameter `parameters::line_distance_threshold`.}
  *   \cgalParamSectionEnd
  *
  *   \cgalParamSectionBegin{do_optimize}
  *     \cgalParamDescription{a boolean indicating if the optimization step is called.}
  *     \cgalParamDefault{`true`}
  *   \cgalParamSectionEnd
  *
  *   \cgalParamSectionBegin{maximum_number_of_iterations}
  *     \cgalParamDescription{the maximum number of iterations of the gradient descent during the optimization step.}
  *     \cgalParamDefault{`Size(20)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{start_step_size}
  *     \cgalParamDescription{the step size at the first iteration of the gradient descent during the optimization step,
  *                           expressed as the size of the maximum voxel edge size.}
  *     \cgalParamDefault{`FT(1.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{end_step_size}
  *     \cgalParamDescription{the step size at the last iteration of the gradient descent during the optimization step,
  *                           expressed as the size of the maximum voxel edge size.}
  *     \cgalParamDefault{`FT(0.125)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{mininmum_energy_delta}
  *     \cgalParamDescription{the minimum energy change to stop the gradient descent iterations during the optimization step.}
  *     \cgalParamDefault{`FT(1.e-3)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{collapse_distance}
  *     \cgalParamDescription{the distance to collapse adjacent points in a line during the gradient descent of the optimization step,
  *                           expressed as the size of the maximum voxel edge size.}
  *     \cgalParamDefault{`FT(0.5)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{smoothing_factor}
  *     \cgalParamDescription{the smoothing factor of the energy of the optimization step.
  *                           0 means no smoothing,
  *                           1 means that the energy will consider smoothing with the same weight
  *                           as the displacement toward the sharp features of the surface.}
  *     \cgalParamDefault{`FT(1.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{normal_refinement_distance}
  *     \cgalParamDescription{the distance to refine the normals of elements near the sharp features during the optimization step,
  *                           expressed as the size of the maximum voxel edge size.}
  *     \cgalParamDefault{`FT(4.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{plane_detection_distance}
  *     \cgalParamDescription{the distance to collect elements near the sharp features
  *                           to determine the adjacent planes during the optimization step,
  *                           expressed as the size of the maximum voxel edge size.}
  *     \cgalParamDefault{`FT(4.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{normal_estimator}
  *     \cgalParamDescription{an estimator to evaluate the normals on elements.}
  *     \cgalParamDefault{`CGAL::Feature_graph::Image_AmbrosioTortorelli(image).normal_estimator()`}
  *   \cgalParamSectionEnd
  *
  *   \cgalParamSectionBegin{point_to_element_output_map}
  *     \cgalParamDescription{an output property map that will be filled if supplied.
  *     It enables the user to retrieve the surface element where a feature graph point is embedded.
  *     It must be a model of `WritablePropertyMap`,
  *     and must implement `put(output_pmap, feature_vertex_descriptor, face_descriptor)`
  *     where the key type is the `vertex_descriptor` of the feature graph, and the value type is
  *     a descriptor of an element of the surface, with type `std::size_t`.
  *     The face descriptor correspond to a face element with the tag `CGAL::DimensionTag<2>`}
  *     \cgalParamDefault{`parameters::default()`}
  *   \cgalParamSectionEnd
  * \cgalNamedParamsEnd
  *
  * \returns a graph model of `VertexAndEdgeListGraph`
  * containing the constructed features.
  */
  template <typename CGAL_NP_TEMPLATE_PARAMETERS>
  unspecified_type operator()(const CGAL::Image_3& image, const CGAL_NP_CLASS& np = parameters::default_values()) const;

  /// @}
};

/*!
* \ingroup PkgFeatureGraphDetector
*
* Detector for sharp feature extraction in surfaces.
* It constructs a feature graph model of `VertexAndEdgeListGraph`.
*
* \tparam K the requirements for the geometric objects, model of `Kernel`.
*
* \sa `Detect_sharp_features_on_labeled_image`
*/
template <typename K>
struct Detect_sharp_features_on_surface
{
  /// \name Types
  /// @{

  /*!
  * Output feature graph point type
  */
  typedef typename K::Point_3 Point_3;

  /*!
  * Natural number type.
  */
  typedef std::size_t Size;

  /*!
  * Numerical type.
  */
  typedef typename K::FT FT;

  /// @}

  /// \name Detector
  /// @{

  /*!
  * \brief detects and constructs the feature graph from a surface
  *
  * The feature graph has polylines that lie on the
  * sharp features of the surface.
  *
  * \tparam PolygonMesh a model of `FaceListGraph` that represents a surface mesh.
  *
  * \param pmesh a polygon mesh from where the sharp features are extracted.
  *
  * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
  * \cgalNamedParamsBegin
  *   \cgalParamSectionBegin{sharpness_estimator}
  *     \cgalParamDescription{an estimator model of `SharpnessEstimator`.
  *       See \ref PkgFeatureGraphSharpnessEstimator for available estimators.}
  *     \cgalParamDefault{`CGAL::Feature_graph::Surface_sharpness_estimator(pmesh)`}
  *   \cgalParamSectionEnd
  *
  *   \cgalParamSectionBegin{line_distance_threshold}
  *     \cgalParamDescription{the threshold on the distance between lines.
  *         In the regularization step, if the maximum distance of a line to another line
  *         is less than this threshold, then the line is collapsed.
  *         It can be a constant or a scalar field that return the threshold at a point on the surface.
  *         If it is a scalar field, it must implement
  *         <UL>
  *           <LI> `template <typename DimensionTag, typename Descriptor>`
  *           <BR> `FT operator()(const Point_3& point_in_space, const Descriptor& element_descriptor, const PolygonMesh& pmesh) const`
  *         </UL>
  *         <LI> The type `DimensionTag` represents the element type, it can either be `CGAL::Dimension_tag<0>`, `CGAL::Dimension_tag<1>` or `CGAL::Dimension_tag<2>`,
  *         which represent either a vertex element, an edge element, or a face element respectively.
  *         <LI> The type `Descriptor` is either a `vertex_descriptor`, `halfedge_descriptor` or `face_descriptor` and
  *         is set respectivelly to the DimensionTag.}
  *     \cgalParamDefault{`FT(0.0)`}
  *     \cgalParamExtra{This parameter sets a baseline distance and is overridden
  *                     by the parameters `parameters::isthmus_line_distance_threshold`, `parameters::simple_line_distance_threshold` and `parameters::corner_line_distance_threshold`
  *                     if they are set.}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{isthmus_line_distance_threshold}
  *     \cgalParamDescription{the threshold on the distance between an isthmus line and any line.
  *         An isthmus line is bounded by a corner that is incident to only the same line.
  *         In the regularization step, if the maximum distance of a line to another line
  *         is less than this threshold, then the line is collapsed.
  *         It can be a constant or a scalar field that return the threshold at a point on the surface.
  *         If it is a scalar field, it must implement
  *         <UL>
  *           <LI> `template <typename DimensionTag, typename Descriptor>`
  *           <BR> `FT operator()(const Point_3& point_in_space, const Descriptor& element_descriptor, const PolygonMesh& pmesh) const`
  *         </UL>
  *         See the description of `parameters::line_distance_threshold` for more details on this scalar field.}
  *     \cgalParamDefault{`FT(-1.0)`}
  *     \cgalParamExtra{Negative values implies the use of the parameter `parameters::line_distance_threshold`.}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{simple_line_distance_threshold}
  *     \cgalParamDescription{the threshold on the distance between a simple line and any line.
  *         A simple line is a line that can be removed without splitting the graph of its neighbors.
  *         In the regularization step, if the maximum distance of a line to another line
  *         is less than this threshold, then the line is collapsed.
  *         It can be a constant or a scalar field that return the threshold at a point on the surface.
  *         If it is a scalar field, it must implement
  *         <UL>
  *           <LI> `template <typename DimensionTag, typename Descriptor>`
  *           <BR> `FT operator()(const Point_3& point_in_space, const Descriptor& element_descriptor, const PolygonMesh& pmesh) const`
  *         </UL>
  *         See the description of `parameters::line_distance_threshold` for more details on this scalar field.}
  *     \cgalParamDefault{`FT(-1.0)`}
  *     \cgalParamExtra{Negative values implies the use of the parameter `parameters::line_distance_threshold`.}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{corner_line_distance_threshold}
  *     \cgalParamDescription{the threshold on the distance between a corner line and any line.
  *         A corner line is neither isthmus nor simple.
  *         In the regularization step, if the maximum distance of a line to another line
  *         is less than this threshold, then the line is collapsed.
  *         It can be a constant or a scalar field that return the threshold at a point on the surface.
  *         If it is a scalar field, it must implement
  *         <UL>
  *           <LI> `template <typename DimensionTag, typename Descriptor>`
  *           <BR> `FT operator()(const Point_3& point_in_space, const Descriptor& element_descriptor, const PolygonMesh& pmesh) const`
  *         </UL>
  *         See the description of `parameters::line_distance_threshold` for more details on this scalar field.}
  *     \cgalParamDefault{`FT(-1.0)`}
  *     \cgalParamExtra{Negative values implies the use of the parameter `parameters::line_distance_threshold`.}
  *   \cgalParamSectionEnd
  *
  *   \cgalParamSectionBegin{do_optimize}
  *     \cgalParamDescription{a boolean indicating if the optimization step is called.}
  *     \cgalParamDefault{`true`}
  *   \cgalParamSectionEnd
  *
  *   \cgalParamSectionBegin{maximum_number_of_iterations}
  *     \cgalParamDescription{the maximum number of iterations of the gradient descent during the optimization step.}
  *     \cgalParamDefault{`Size(20)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{start_step_size}
  *     \cgalParamDescription{the step size at the first iteration of the gradient descent during the optimization step.}
  *     \cgalParamDefault{`FT(0.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{end_step_size}
  *     \cgalParamDescription{the step size at the last iteration of the gradient descent during the optimization step.}
  *     \cgalParamDefault{`FT(0.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{mininmum_energy_delta}
  *     \cgalParamDescription{the minimum energy change to stop the gradient descent iterations during the optimization step.}
  *     \cgalParamDefault{`FT(1.e-3)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{collapse_distance}
  *     \cgalParamDescription{the distance to collapse adjacent points in a line during the gradient descent of the optimization step.}
  *     \cgalParamDefault{`FT(0.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{smoothing_factor}
  *     \cgalParamDescription{the smoothing factor of the energy of the optimization step.
  *                           0 means no smoothing,
  *                           1 means that the energy will consider smoothing with the same weight
  *                           as the displacement toward the sharp features of the surface.}
  *     \cgalParamDefault{`FT(1.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{normal_refinement_distance}
  *     \cgalParamDescription{the distance to refine the normals of elements near the sharp features during the optimization step.}
  *     \cgalParamDefault{`FT(0.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{plane_detection_distance}
  *     \cgalParamDescription{the distance to collect elements near the sharp features
  *                           to determine the adjacent planes during the optimization step.}
  *     \cgalParamDefault{`FT(0.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{normal_estimator}
  *     \cgalParamDescription{an estimator to evaluate the normals on elements.}
  *     \cgalParamDefault{`CGAL::Feature_graph::Surface_normal_estimator(pmesh)`}
  *   \cgalParamSectionEnd
  *
  *   \cgalParamSectionBegin{point_to_element_output_map}
  *     \cgalParamDescription{an output property map that will be filled if supplied.
  *     It enables the user to retrieve the surface element where a feature graph point is embedded.
  *     It must be a model of `WritablePropertyMap`,
  *     and must implement `put(output_pmap, feature_vertex_descriptor, face_descriptor)`
  *     where the key type is the `vertex_descriptor` of the feature graph, and the value type is
  *     the `facet_descriptor` of the polygon mesh.}
  *     \cgalParamDefault{`parameters::default()`}
  *   \cgalParamSectionEnd
  * \cgalNamedParamsEnd
  *
  * \returns a graph model of `VertexAndEdgeListGraph`
  * containing the constructed features.
  */
  template <typename PolygonMesh, typename CGAL_NP_TEMPLATE_PARAMETERS>
  unspecified_type operator()(const PolygonMesh& pmesh, const CGAL_NP_CLASS& np = parameters::default_values()) const;

  /// @}
};

} /* namespace Feature_graph */

} /* namespace CGAL */