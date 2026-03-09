namespace CGAL {


/*!
\ingroup PkgMesh3Domains

\brief The class `Labeled_mesh_domain_3` implements indexed domains.

This class is a model of concept `MeshDomain_3`.

Any boundary facet is labeled <a,b>, with a<b, where a and b are the
tags of its incident subdomains.
Thus, a boundary facet of the domain is labeled <0,b>, where b!=0.

This class includes a \link Labeling_function `Labeling_function` \endlink that provides the index of the subdomain in which any
query point lies. An intersection between a segment and bounding
surfaces is detected when both segment endpoints are associated with different
values of subdomain indices. The intersection is then constructed by bisection.
The bisection stops when the query segment is shorter than an error bound
`e` given by the product of the
length of the diagonal of the bounding box (in world coordinates), or the radius of the bounding sphere, and
a relative error bound passed as argument to the constructor of `Labeled_mesh_domain_3`.

This class has a constructor taking a labeling function. It has also three
static template member functions that act as named constructors:
<ul>
<li>`create_gray_image_mesh_domain()`, to create a domain from a 3D gray image,
<li>`create_labeled_image_mesh_domain()`, to create a domain from a 3D labeled image, and
<li>`create_implicit_mesh_domain()`, to create a domain from an implicit function.
</ul>

\tparam BGT is a geometric traits class that provides
the basic operations to implement
intersection tests and intersection computations
through a bisection method. This parameter must be instantiated
with a model of the concept `BisectionGeometricTraits_3`.

\tparam SubdomainIndex is the type of the indices of the subdomains.
\tparam SurfacePatchIndex is the type of the indices of the surface patches.


\cgalModels{MeshDomain_3}

\sa `CGAL::Implicit_multi_domain_to_labeling_function_wrapper`
\sa `CGAL::make_mesh_3()`

*/
template<class BGT,
         class SubdomainIndex = int,
         class SurfacePatchIndex = std::pair<SubdomainIndex,
                                             SubdomainIndex> >
class Labeled_mesh_domain_3
{
public:
  //-------------------------------------------------------
  // Index Types
  //-------------------------------------------------------
  // Type of indexes for cells of the input complex

/// \name Types
///@{



  /// The subdomain index of this model of `MeshDomain_3`
  typedef Subdomain_index_                  Subdomain_index;



  /// \brief The type of object that stores the function using type-erasure.
  /// \details A labeling function `f` must return `0` if the point is not located in any subdomain. The return type of labeling functions is an integer.
  ///
  /// Let `p` be a point.
  /// <ul>
  /// <li>`f(p)=0` means that `p` is outside the domain.</li>
  /// <li>`f(p)=a`, `a!=0` means that `p` is inside subdomain `a`.</li>
  /// </ul>
  ///
  /// The function type must be a model of the concept `Callable` with signature
  /// `Subdomain_index(const %Point_3&)`: it can be a function, a function object, a lambda expression...
  /// that takes a `%Point_3` as argument, and returns a type convertible to `Subdomain_index`.

  typedef std::function< Subdomain_index(const Point_3 &)> Labeling_function;
///@}

/// \name Types imported from the geometric traits class
///@{
  /// The point type of the geometric traits class
  typedef typename BGT::Point_3      Point_3;
  /// The sphere type of the geometric traits class
  typedef typename BGT::Sphere_3     Sphere_3;
  /// The iso-cuboid type of the geometric traits class
  typedef typename BGT::Iso_cuboid_3 Iso_cuboid_3;
  /// The number type (a field type) of the geometric traits class
  typedef typename BGT::FT           FT;
///@}


/// \name Creation
/// @{
  /*!  \brief Construction from a function, a bounding object and a relative error bound.
   *
   * \tparam Function a type compatible with \link Labeling_function `Labeling_function`\endlink.
   *   The class `CGAL::Implicit_multi_domain_to_labeling_function_wrapper` is a good candidate for this template parameter
   *   if there are several components to mesh.
   * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   * \tparam BoundingObject either a bounding sphere (of type `Sphere_3`), a bounding box (type `Bbox_3` or `Iso_cuboid_3`)
   *
   * \param function the labeling function
   * \param bounding_object the bounding object bounding the meshable space.
   * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{relative_error_bound}
   *      \cgalParamDescription{the relative error bound used to compute intersection points between the implicit surface and query segments.
   *                            The bisection is stopped when the length of the intersected segment is less than the product
   *                            of `relative_error_bound` by the diameter of the bounding object.}
   *      \cgalParamDefault{FT(1e-3)}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   *
   * \cgalHeading{Example}
   * From the example (\ref Mesh_3/mesh_implicit_domains_2.cpp):
   * \snippet Mesh_3/mesh_implicit_domains_2.cpp Domain creation
   */
  template<typename Function, typename Bounding_object, typename CGAL_NP_TEMPLATE_PARAMETERS>
  Labeled_mesh_domain_3(const Function& function,
                        const Bounding_object& bounding_object,
                        const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                        , typename std::enable_if<!is_named_function_parameter<Function>>::type* = nullptr
#endif // DOXYGEN_RUNNING
                        );
///@}

/// \name Creation of domains from 3D images
/// @{
  /*!
   * \brief Construction from a 3D gray image
   *
   * This static method is a <em>named constructor</em>. It constructs a domain
   * described by a 3D gray image. A 3D gray image is a grid of voxels,
   * where each voxel is associated with a gray level value.  Unless otherwise specified by the parameter `image_values_to_subdom_indices`, the domain to
   * be discretized is the union of voxels that lie inside a surface
   * described by an isolevel value, called \a isovalue. The voxels lying
   * inside the domain have gray level values that are larger than the
   * isovalue.
   *
   * The value of voxels is interpolated to a gray level value at any query point.
   *
   * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   *
   * \param image_ the input 3D image.
   * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{iso_value}
   *     \cgalParamDescription{the isovalue, inside
   *                           `image`, of the surface describing the boundary of the object to be
   *                            meshed.}
   *     \cgalParamDefault{0}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{image_values_to_subdomain_indices}
   *     \cgalParamDescription{a function or a function object, compatible with the signature
   *                           `Subdomain_index(double)`. This function returns the subdomain index
   *                           corresponding to a pixel value. If this parameter is used, then the
   *                           parameter `iso_value` is ignored.}
   *     \cgalParamDefault{Null_functor()}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{value_outside}
   *      \cgalParamDescription{the value attached to voxels
   *                            outside of the domain to be meshed. It should be lower than
   *                           `iso_value`.}
   *      \cgalParamDefault{0}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{relative_error_bound}
   *      \cgalParamDescription{ is the relative error
   *                             bound, relative to the diameter of the box of the image.}
   *      \cgalParamDefault{FT(1e-3)}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   *
   * \cgalHeading{Examples}
   *
   * From the example (\ref Mesh_3/mesh_3D_gray_image.cpp):
   *
   * \snippet Mesh_3/mesh_3D_gray_image.cpp Domain creation
   *
   * From the example (\ref Mesh_3/mesh_3D_gray_vtk_image.cpp):
   *
   * \snippet Mesh_3/mesh_3D_gray_vtk_image.cpp Domain creation
   *
   */
  template<typename CGAL_NP_TEMPLATE_PARAMETERS>
  static Labeled_mesh_domain_3 create_gray_image_mesh_domain(const CGAL::Image_3& image_, const CGAL_NP_CLASS& np = parameters::default_values());

  /*!
   * \brief Construction from a 3D labeled image
   *
   * This static method is a <em>named constructor</em>. It constructs a
   * domain described by a 3D labeled image. A 3D labeled image is a grid
   * of voxels, where each voxel is associated with an index (a subdomain
   * index) characterizing the subdomain in which the voxel lies. The
   * domain to be discretized is the union of voxels that have non-zero
   * values.
   *
   * \returns either a `Labeled_mesh_domain_3`,
   *   or a `Mesh_domain_with_polyline_features_3<Labeled_mesh_domain_3>`
   *   depending on whether one or more of the named parameters
   *   `features_detector` and `input_features` are provided.
   *
   * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   * \param image_ the input 3D image.
   * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{weights}
   *     \cgalParamDescription{a reference to an input 3D image that provides
   *                           weights associated to each voxel (the word type is `unsigned char`,
   *                           and the voxels values are integers between 0 and 255).
   *                           The weights image can be generated with `CGAL::Mesh_3::generate_label_weights()`.
   *                           Its dimensions must be the same as the dimensions of `parameters::image`.}
   *     \cgalParamDefault{CGAL::Image_3()}
   *     \cgalParamType{CGAL::Image_3&}
   *     \cgalParamExtra{if `features_detector` is provided, `weights` should be modified accordingly.
   *                     The available functors described in See \ref PkgMesh3FeatureDetection
   *                     implement the necessary modifications.}
   *     \cgalParamExtra{if `input_features` is provided, `weights` should be modified accordingly
   *                     to keep consistency of the output `MeshDomainWithFeatures_3`}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{value_outside}
   *     \cgalParamDescription{the value attached to voxels
   *                           outside of the domain to be meshed. It should be lower than
   *                           `iso_value`.}
   *     \cgalParamDefault{0}
   *   \cgalParamNEnd
   *
   *   \cgalParamNBegin{relative_error_bound}
   *     \cgalParamDescription{ is the relative error
   *                            bound, relative to the diameter of the box of the image.}
   *     \cgalParamDefault{FT(1e-3)}
   *   \cgalParamNEnd
   *
   *   \cgalParamNBegin{features_detector}
   *    \cgalParamDescription{ a functor that implements
   *      `std::vector<std::vector<Point>> operator()(const Image_3& img) const`,
   *      and `std::vector<std::vector<Point>> operator()(const Image_3& img, Image_3& weights) const`,
   *      where `%Point` matches the mesh domain point type,
   *      that both return a range of detected polyline features for feature protection.
   *      Only one implementation is used, depending on whether the named parameter `weights`
   *      is provided or not.
   *      Polyline features are added to the domain for further feature protection.
   *      See \ref PkgMesh3FeatureDetection for available functors.}
   *    \cgalParamDefault{CGAL::Null_functor()}
   *    \cgalParamExtra{The return type of the function depends on whether this parameter
   *                    or `input_features` are provided or not.}
   *    \cgalParamExtra{If `weights` is provided, it must either be adapted to the detected features,
   *                    or postprocessed during feature detection to keep consistency
   *                    of the output `MeshDomainWithFeatures_3`.
   *                    Available functors implement the necessary modifications.}
   *   \cgalParamNEnd
   *
   *   \cgalParamNBegin{input_features}
   *    \cgalParamDescription{ a `Range` of polyline features, represented as `Range`s of `Point_3`.
   *         Polyline features are added to the domain for further feature protection.
   *         Input polyline features must be different from the detected features
   *         and can intersect only at vertices, if they do. Otherwise,
   *         the meshing process may not terminate.}
   *    \cgalParamDefault{`std::vector<std::vector<Point_3>>()`}
   *    \cgalParamExtra{The return type of the function depends on whether this parameter
                        or `input_features` are provided or not.}
   *    \cgalParamExtra{It is recommended to pass a const-reference for this parameter,
   *                    possibly using `std::cref(polylines_range)` to avoid useless copies.}
   *    \cgalParamExtra{If `weights` is provided, it must be adapted to the input features,
   *                    to keep consistency of the output `MeshDomainWithFeatures_3`}
   *   \cgalParamNEnd
   *
   * \cgalNamedParamsEnd
   *
   * \cgalHeading{Example}
   *
   * From the example (\ref Mesh_3/mesh_3D_image.cpp):
   *
   * \snippet Mesh_3/mesh_3D_image.cpp Domain creation
   *
   * From the example (\ref Mesh_3/mesh_3D_weighted_image.cpp),
   * where the labeled image is used with a precomputed 3D image of weights :
   *
   * \snippet Mesh_3/mesh_3D_weighted_image.cpp Domain creation
   *
   * From the example (\ref Mesh_3/mesh_3D_image_with_detection_of_features.cpp)
   * where the features are detected in `image`:
   *
   * \snippet Mesh_3/mesh_3D_image_with_detection_of_features.cpp Domain creation
   *
   * From the example (\ref Mesh_3/mesh_3D_image_with_input_features.cpp)
   * where the features are provided by the user:
   *
   * \snippet Mesh_3/mesh_3D_image_with_input_features.cpp Domain creation
   *
   */
  template<typename CGAL_NP_TEMPLATE_PARAMETERS>
  static auto
  create_labeled_image_mesh_domain(const CGAL::Image_3& image_, const CGAL_NP_CLASS& np = parameters::default_values());
/// @}


/// \name Creation of domains from implicit functions
/// @{

  /*!
   * \brief Construction from an implicit function
   *
   * This static method is a <em>named constructor</em>. It constructs a domain
   * whose bounding surface is described implicitly as the zero level set of a
   * function.  The domain to be discretized is assumed to be the domain where
   * the function has negative values.
   *
   * \tparam Function a type compatible with the signature `FT(Point_3)`: it takes a point as argument,
   *                  and returns a scalar value. That object must be model of `CopyConstructible`
   * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   * \tparam BoundingObject either a bounding sphere (of type `Sphere_3`), a bounding box
   *         (type `Bbox_3` or `Iso_cuboid_3`) which is required to circumscribe
   *                         the surface and to have its center inside the domain.
   *
   * \param function the implicit function
   * \param bounding_object object bounding the meshable domain and its center is inside the domain.
   * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{relative_error_bound}
   *     \cgalParamDescription{ is the relative error
   *                            bound, relative to the diameter of the box of the image.}
   *     \cgalParamDefault{FT(1e-3)}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   *
   * \cgalHeading{Examples}
   *
   * From the example (\ref Mesh_3/mesh_implicit_sphere.cpp):
   *
   * \snippet Mesh_3/mesh_implicit_sphere.cpp Domain creation
   *
   * From the example (\ref Mesh_3/mesh_implicit_sphere_variable_size.cpp):
   *
   * \snippet Mesh_3/mesh_implicit_sphere_variable_size.cpp Domain creation
   *
   */
  template<typename Function, typename BoundingObject, typename CGAL_NP_TEMPLATE_PARAMETERS>
  static Labeled_mesh_domain_3 create_implicit_mesh_domain(const Function& function,
                                                           const BoundingObject& bounding_object,
                                                           const CGAL_NP_CLASS& np = parameters::default_values());
/// @}
  /*
   * Constructs  a set of `n` points on the surface, and output them to
   *  the output iterator `pts` whose value type is required to be
   *  `std::pair<Points_3, Index>`.
   */
  struct Construct_initial_points
  {
    Construct_initial_points(const Labeled_mesh_domain_3& domain);

    template<class OutputIterator>
    OutputIterator operator()(OutputIterator pts, const int n = 12) const;
  };

  // Returns Construct_initial_points object
  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }

  /*
   * Returns a bounding box of the domain
   */
  Bbox_3 bbox() const {
    return this->bbox_.bbox();
  }

  /*
   * Returns `true` if point `p` is in the domain. If `p` is in the
   *  domain, the parameter index is set to the index of the subdomain
   *  including `p`. It is set to the default value otherwise.
   */
  struct Is_in_domain
  {
    Is_in_domain(const Labeled_mesh_domain_3& domain);

    Subdomain operator()(const Point_3& p) const;
  };

  // Returns Is_in_domain object
  Is_in_domain is_in_domain_object() const;

  /*
   * Returns `true` if the element `type` intersect properly any of the
   * surface patches describing the either the domain boundary or some
   * subdomain boundary.
   * `Type` is either `Segment_3`, `Ray_3` or `Line_3`.
   * Parameter index is set to the index of the intersected surface patch
   * if `true` is returned and to the default `Surface_patch_index`
   * value otherwise.
   */
  struct Do_intersect_surface
  {
    Do_intersect_surface(const Labeled_mesh_domain_3& domain);

    Surface_patch operator()(const Segment_3& s) const;
    Surface_patch operator()(const Ray_3& r) const;
    Surface_patch operator()(const Line_3& l) const;
  };

  // Returns Do_intersect_surface object
  Do_intersect_surface do_intersect_surface_object() const;

  /*
   * Returns a point in the intersection of the primitive `type`
   * with some boundary surface.
   * `Type` is either `Segment_3`, `Ray_3` or `Line_3`.
   */
  struct Construct_intersection
  {
    Construct_intersection(const Labeled_mesh_domain_3& domain);
    Intersection operator()(const Segment_3& s) const;
    Intersection operator()(const Ray_3& r) const;
    Intersection operator()(const Line_3& l) const;
  };

  // Returns Construct_intersection object
  Construct_intersection construct_intersection_object() const;

  /*
   * Returns the index to be stored in a vertex lying on the surface identified
   * by `index`.
   */
  Index index_from_surface_patch_index(const Surface_patch_index& index) const;

  /*
   * Returns the index to be stored in a vertex lying in the subdomain
   * identified by `index`.
   */
  Index index_from_subdomain_index(const Subdomain_index& index) const;

  /*
   * Returns the `Surface_patch_index` of the surface patch
   * where lies a vertex with dimension 2 and index `index`.
   */
  Surface_patch_index surface_patch_index(const Index& index) const;

  /*
   * Returns the index of the subdomain containing a vertex
   *  with dimension 3 and index `index`.
   */
  Subdomain_index subdomain_index(const Index& index) const;
};  // end class Labeled_mesh_domain_3

}  // end namespace CGAL
