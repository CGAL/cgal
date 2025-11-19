// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sven Oesau

#ifndef CGAL_POISSON_MESH_DOMAIN_3_H
#define CGAL_POISSON_MESH_DOMAIN_3_H

#include <CGAL/license/Poisson_surface_reconstruction_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Poisson_reconstruction_function.h>

namespace CGAL {

/*!
\ingroup PkgPoissonSurfaceReconstruction3Ref

\brief The class `Poisson_mesh_domain_3` derives from `Labeled_mesh_domain_3` for the handling of `Poisson_reconstruction_function`.

 This class has a constructor taking a labeling function. It has also a static template member
 function that acts as named constructor:
 <ul><li>`create_Poisson_mesh_domain()`</li>, to create a domain from a `Poisson_reconstruction_function`</ul>

\tparam BGT is a geometric traits class that provides
    the basic operations to implement intersection tests and intersection computations through a bisection
    method. This parameter must be instantiated with a model of the concept `BisectionGeometricTraits_3`.

\cgalModels{MeshDomain_3}

\sa `CGAL::Labeled_mesh_domain_3`
\sa `CGAL::make_mesh_3()`
*/
template<class BGT>
class Poisson_mesh_domain_3
#ifndef DOXYGEN_RUNNING
    : public Labeled_mesh_domain_3<BGT>
#endif
{
public:
  using Base = Labeled_mesh_domain_3<BGT>;
  typedef typename Base::Subdomain Subdomain;
  typedef typename Base::Subdomain_index Subdomain_index;
  typedef typename Base::Surface_patch_index Surface_patch_index;
  typedef typename Base::Intersection Intersection;

  // Type of indexes for cells of the input complex
  typedef std::optional<Surface_patch_index> Surface_patch;

  // Type of indexes to characterize the lowest dimensional face of the input
  // complex on which a vertex lie
  typedef typename CGAL::Mesh_3::internal::Index_generator<Subdomain_index, Surface_patch_index>::Index Index;

  // Geometric object types
#ifdef DOXYGEN_RUNNING
/// \name Types imported from the geometric traits class
  ///@{
  /// The point type of the geometric traits class
  typedef typename Geom_traits::Point_3 Point_3;
  /// The sphere type of the geometric traits class
  typedef typename Geom_traits::Sphere_3 Sphere_3;
  /// The iso-cuboid type of the geometric traits class
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid_3;
  /// The bounding box type
  typedef CGAL::Bbox_3 Bbox_3;
  /// The number type (a field type) of the geometric traits class
  typedef typename Geom_traits::FT FT;
  /// The ray type of the geometric traits class
  typedef typename Geom_traits::Ray_3 Ray_3;
  /// The line type of the geometric traits class
  typedef typename Geom_traits::Line_3 Line_3;
  /// The segment type of the geometric traits class
  typedef typename Geom_traits::Segment_3 Segment_3;
  /// The Poisson function type
  typedef CGAL::Poisson_reconstruction_function<Geom_traits> Function;
  ///@}
#else
  /// The point type of the geometric traits class
  typedef typename BGT::Point_3 Point_3;
  /// The sphere type of the geometric traits class
  typedef typename BGT::Sphere_3 Sphere_3;
  /// The iso-cuboid type of the geometric traits class
  typedef typename BGT::Iso_cuboid_3 Iso_cuboid_3;
  /// The bounding box type
  typedef CGAL::Bbox_3 Bbox_3;
  /// The number type (a field type) of the geometric traits class
  typedef typename BGT::FT FT;
  /// The ray type of the geometric traits class
  typedef typename BGT::Ray_3 Ray_3;
  /// The line type of the geometric traits class
  typedef typename BGT::Line_3 Line_3;
  /// The segment type of the geometric traits class
  typedef typename BGT::Segment_3 Segment_3;
  /// The Poisson function type
  typedef CGAL::Poisson_reconstruction_function<BGT> Function;
#endif

  Function poisson_function;

/// \name Creation
/// @{
  /*!  \brief Construction from a function, a bounding object and a relative error bound.
   *
   * \tparam Bounding_object either a bounding sphere (of type `Sphere_3`), a bounding box (type `Bbox_3`),
   *                         or a bounding `Iso_cuboid_3`
   * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   *
   * \param function the Poisson reconstruction function
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
   */
  template<typename Bounding_object, typename CGAL_NP_TEMPLATE_PARAMETERS>
  Poisson_mesh_domain_3(const Function& function,
    const Bounding_object& bounding_object,
    const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
    , typename std::enable_if<!is_named_function_parameter<Function>>::type* = nullptr
#endif // DOXYGEN_RUNNING
  )
    : Base(make_implicit_to_labeling_function_wrapper<BGT>(function), bounding_object, np),
      poisson_function(function)
  {}

  /*!  \brief Construction from a function, a bounding object and a relative error bound.
   *
   * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   *
   * \param function the Poisson reconstruction function
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
   */
  template<typename CGAL_NP_TEMPLATE_PARAMETERS>
  Poisson_mesh_domain_3(const Function & function,
    const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
    , typename std::enable_if<!is_named_function_parameter<Function>>::type * = nullptr
#endif // DOXYGEN_RUNNING
  )
    : Base(make_implicit_to_labeling_function_wrapper<BGT>(function), function.bounding_sphere(), np),
      poisson_function(function)
  {}
///@}

#ifndef DOXYGEN_RUNNING
  template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
  Poisson_mesh_domain_3(const CGAL_NP_CLASS& np)
    : Base(np)
  {}

  // Overload handling parameters passed with operator=

  template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
           typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
           typename ... NP>
  Poisson_mesh_domain_3(const Function& function,
                        const CGAL_NP_CLASS_1& np1,
                        const CGAL_NP_CLASS_2& np2,
                        const NP& ... nps)
    : Base(internal_np::combine_named_parameters(
          CGAL::parameters::function(make_implicit_to_labeling_function_wrapper<BGT>(function)), np1, np2, nps...)),
      poisson_function(function)
  {}
#endif

/// \name Creation of domains from Poisson implicit functions
/// @{
  /*! \brief Construction from a Poisson implicit function
   *
   * This static method is a <em>named constructor</em>. It constructs a domain
   * whose bounding surface is described implicitly as the zero level set of a
   * function.  The domain to be discretized is assumed to be the domain where
   * the function has negative values.
   *
   * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   * \tparam Bounding_object either a bounding sphere (of type `Sphere_3`), a bounding box (type `Bbox_3`),
   *                         or a bounding `Iso_cuboid_3` which is required to circumscribe
   *                         the surface and to have its center inside the domain.
   *
   * \param function the Poisson reconstruction function
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
   */
  template<typename Bounding_object, typename CGAL_NP_TEMPLATE_PARAMETERS>
  static Poisson_mesh_domain_3 create_Poisson_mesh_domain(const Function& function,
                                                          const Bounding_object& bounding_object,
                                                          const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                                                          ,typename std::enable_if<!is_named_function_parameter<Function>>::type* = nullptr
#endif
)
  {
    using parameters::get_parameter;
    using parameters::choose_parameter;

    FT relative_error_bound_ = choose_parameter(get_parameter(np, internal_np::error_bound), FT(1e-3));
    CGAL::Random* p_rng_ = choose_parameter(get_parameter(np, internal_np::rng), nullptr);
    auto null_subdomain_index_ = choose_parameter(get_parameter(np, internal_np::null_subdomain_index_param), Null_functor());
    auto construct_surface_patch_index_ = choose_parameter(get_parameter(np, internal_np::surface_patch_index), Null_functor());

    return Poisson_mesh_domain_3(function,
             bounding_object,
             CGAL::parameters::relative_error_bound(relative_error_bound_)
             .function(make_implicit_to_labeling_function_wrapper<BGT>(function))
             .p_rng(p_rng_)
             .null_subdomain_index(Base::create_null_subdomain_index(null_subdomain_index_))
             .construct_surface_patch_index(Base::create_construct_surface_patch_index(construct_surface_patch_index_)));
  }
/// @}
#ifndef DOXYGEN_RUNNING
  template<typename CGAL_NP_TEMPLATE_PARAMETERS>
  static Poisson_mesh_domain_3 create_Poisson_mesh_domain(const CGAL_NP_CLASS& np) {
    using parameters::get_parameter;

    static_assert(!parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::function_param_t>::value, "Value for required parameter not found");
    static_assert(!parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::bounding_object_param_t>::value, "Value for required parameter not found");

    return create_Poisson_mesh_domain(get_parameter(np, internal_np::function_param),
                                       get_parameter(np, internal_np::bounding_object_param),
                                       np);
  }

  // Overload handling parameters passed with operator=
  template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
           typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
           typename ... NP>
  static Poisson_mesh_domain_3 create_Poisson_mesh_domain(const CGAL_NP_CLASS_1& np1,
                                                          const CGAL_NP_CLASS_2& np2,
                                                          const NP& ... nps)
  {
    return create_Poisson_mesh_domain(internal_np::combine_named_parameters(np1, np2, nps...));
  }
#endif

  /*
   * Returns a point in the intersection of the primitive `type`
   * with some boundary surface.
   * `Type` is either `Segment_3`, `Ray_3` or `Line_3`.
   */
  struct Construct_intersection
  {
    Construct_intersection(const Poisson_mesh_domain_3& domain) : domain_(domain) {}

    Intersection operator()(const Segment_3& s) const
    {
#ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
      CGAL_precondition(r_domain_.do_intersect_surface_object()(s) != std::nullopt);
#endif // NOT CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
      return this->operator()(s.source(), s.target());
    }

    Intersection operator()(const Ray_3& r) const
    {
      return clip_to_segment(r);
    }

    Intersection operator()(const Line_3& l) const
    {
      return clip_to_segment(l);
    }

  private:
    /*
     * Returns a point in the intersection of `[a,b]` with the surface
     *  `a` must be the source point, and `b` the out point. It is important
     * because it drives bisection cuts.
     * Indeed, the returned point is the first intersection of `[a,b]`
     * with a subdomain surface.
     *
     * The difference from the Labeled_mesh_domain_3::Construct_intersection is that
     * the underlying Delaunay triangulation in the Poisson function is used for the bisection.
     */
    Intersection operator()(const Point_3& a, const Point_3& b) const
    {
      // Functors
      typename BGT::Compute_squared_distance_3 squared_distance =
                                      BGT().compute_squared_distance_3_object();
      typename BGT::Construct_midpoint_3 midpoint =
                                      BGT().construct_midpoint_3_object();

      // Non const points
      Point_3 p1 = a;
      Point_3 p2 = b;
      Point_3 mid = midpoint(p1, p2);

      FT value_at_p1, value_at_p2;
      typename Function::Cell_handle c1, c2;
      bool c1_is_inf, c2_is_inf;

      std::tie(value_at_p1, c1, c1_is_inf) = domain_.poisson_function.special_func(p1);
      std::tie(value_at_p2, c2, c2_is_inf) = domain_.poisson_function.special_func(p2);

      Subdomain_index label_at_p1 = (value_at_p1 < 0) ? 1 : 0;
      Subdomain_index label_at_p2 = (value_at_p2 < 0) ? 1 : 0;

      // If both extremities are in the same subdomain,
      // there is no intersection.
      // Should only be able to happen during initial point generation.
      if(label_at_p1 == label_at_p2)
        return Intersection();

      // Else lets find a point (by bisection)
      // Bisection ends when the point is nearer from surface than the error bound
      while(true) {
        if(c1 == c2) {
          if(c1_is_inf) {
            std::cout << "Intersection(): c1 == c2 and inf!" << std::endl;
            return Intersection();
          } else {
            const Surface_patch_index sp_index = domain_.make_surface_index(label_at_p1, label_at_p2);
            const Index index = domain_.index_from_surface_patch_index(sp_index);
            return Intersection(Point_3(ORIGIN + ((value_at_p2 * (p1 - ORIGIN)) - (value_at_p1 * (p2 - ORIGIN))) /
                                                  (value_at_p2 - value_at_p1)), index, 2);
          }
        }
        mid = midpoint(p1, p2);
        // If the two points are enough close, then we return midpoint
        if ( squared_distance(p1, p2) < domain_.squared_error_bound_ )
        {
          CGAL_assertion(value_at_p1 * value_at_p2 <= 0);
          const Surface_patch_index sp_index = domain_.make_surface_index(label_at_p1, label_at_p2);
          const Index index = domain_.index_from_surface_patch_index(sp_index);
          return Intersection(mid, index, 2);
        }

        // Cannot be const: those values are modified below.
        FT value_at_mid;
        typename Function::Cell_handle c_at_mid;
        bool c_is_inf;
        std::tie(value_at_mid, c_at_mid, c_is_inf) = domain_.poisson_function.special_func(mid);
        Subdomain_index label_at_mid = (value_at_mid < 0) ? 1 : 0;

        // Else we must go on
        // Here we consider that p1(a) is the source point. Thus, we keep p1 and
        // change p2 if f(p1)!=f(p2).
        // That allows us to find the first intersection from a of [a,b] with
        // a surface.
        if(label_at_p1 != label_at_mid && !(domain_.null(label_at_p1) && domain_.null(label_at_mid)))
        {
          p2 = mid;
          value_at_p2 = value_at_mid;
          label_at_p2 = label_at_mid;
        }
        else
        {
          p1 = mid;
          value_at_p1 = value_at_mid;
          label_at_p1 = label_at_mid;
        }
      }
    }

    // Clips  `query` to a segment `s`, and call `operator()(s)`
    template<typename Query>
    Intersection clip_to_segment(const Query& query) const
    {
      const auto clipped = CGAL::intersection(query, domain_.bbox_);
      if (clipped)
        if (const Segment_3* s = std::get_if<Segment_3>(&*clipped))
          return this->operator()(*s);

      return Intersection();
    }

    const Poisson_mesh_domain_3& domain_;
  };

  // Returns Construct_intersection object
  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection(*this);
  }

};  // end class Poisson_mesh_domain_3

} // end namespace CGAL

#endif // CGAL_LABELED_MESH_DOMAIN_3_H
