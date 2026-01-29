// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Stéphane Tayeb, Aymeric PELLE


#ifndef CGAL_LABELED_MESH_DOMAIN_3_IMAGE_H
#define CGAL_LABELED_MESH_DOMAIN_3_IMAGE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_3/config.h>

// support for `CGAL::Image_3`
#include <CGAL/Image_3.h>
#include <CGAL/Mesh_3/Image_to_labeled_function_wrapper.h>
#include <CGAL/Mesh_3/Image_plus_weights_to_labeled_function_wrapper.h>
#include <CGAL/Mesh_3/polylines_to_protect_in_image.h>

namespace CGAL {
namespace Mesh_3 {
namespace internal {

// Returns a box enclosing image `im`
inline Bbox_3 compute_bounding_box(const Image_3& im)
{
  return Bbox_3(-1+im.tx(),-1+im.ty(),-1+im.tz(),
                double(im.xdim())*im.vx()+im.tx()+1,
                double(im.ydim())*im.vy()+im.ty()+1,
                double(im.zdim())*im.vz()+im.tz()+1);
}

// Detect_features_in_domain
template<typename Point, typename DetectFunctor>
struct Detect_features_in_domain {

  std::vector<std::vector<Point>>
  operator()(const CGAL::Image_3& image, CGAL::Image_3& weights, DetectFunctor functor) const {
#if defined(BOOST_MSVC) && (BOOST_MSVC < 1910) //before msvc2017
    if(weights.is_valid())
      return functor.operator()<Point>(image, weights);
    else
      return functor.operator()<Point>(image);
#else
    if(weights.is_valid())
      return functor.template operator()<Point>(image, weights);
    else
      return functor.template operator()<Point>(image);
#endif
  }
};
// specialization for `Null_functor`: create the default functor
template<typename Point>
struct Detect_features_in_domain<Point, Null_functor> {
  std::vector<std::vector<Point>>
  operator()(const CGAL::Image_3&, CGAL::Image_3&, Null_functor) const {
    return std::vector<std::vector<Point>>();
  }
};

template<typename Point, typename DetectFunctor>
std::vector<std::vector<Point>>
  detect_features(const CGAL::Image_3& image, CGAL::Image_3& weights, DetectFunctor functor)
{
  Detect_features_in_domain<Point, DetectFunctor> detector;
  return detector(image, weights, functor);
}

template<bool WithFeatures>
struct Add_features_in_domain {
  template<typename MeshDomain, typename InputFeatureRange, typename DetectFunctor>
  void operator()(const CGAL::Image_3&, CGAL::Image_3&, MeshDomain&, const InputFeatureRange&, DetectFunctor)
  {}
};

template<>
struct Add_features_in_domain<true>
{
  template<typename MeshDomain, typename InputFeatureRange, typename DetectFunctor>
  void operator()(const CGAL::Image_3& image,
                  CGAL::Image_3& weights,
                  MeshDomain& domain,
                  const InputFeatureRange& input_features,
                  DetectFunctor functor)
  {
    using P = typename MeshDomain::Point_3;
    auto detected_feature_range
      = CGAL::Mesh_3::internal::detect_features<P>(image, weights, functor);

    CGAL::merge_and_snap_polylines(image, detected_feature_range, input_features);

    if (!input_features.empty())
      domain.add_features(input_features.begin(), input_features.end());
    domain.add_features(detected_feature_range.begin(), detected_feature_range.end());
  }
};

} } // end of Mesh_3::internal


template<class BGT, class Subdomain_index, class Surface_patch_index>
template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_gray_image_mesh_domain(const CGAL::Image_3& image_, const CGAL_NP_CLASS& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;
  auto iso_value_ = choose_parameter(get_parameter(np, internal_np::iso_value_param), 0);
  auto value_outside_ = choose_parameter(get_parameter(np, internal_np::voxel_value), 0);
  FT relative_error_bound_ = choose_parameter(get_parameter(np, internal_np::error_bound), FT(1e-3));
  auto image_values_to_subdomain_indices_ = choose_parameter(get_parameter(np, internal_np::image_subdomain_index), Null_functor());
  CGAL::Random* p_rng_ = choose_parameter(get_parameter(np, internal_np::rng), nullptr);
  auto null_subdomain_index_ = choose_parameter(get_parameter(np, internal_np::null_subdomain_index_param), Null_functor());
  auto construct_surface_patch_index_ = choose_parameter(get_parameter(np, internal_np::surface_patch_index), Null_functor());
  namespace p = CGAL::parameters;
  return Labeled_mesh_domain_3
            (p::function = create_gray_image_wrapper
                     (image_,
                      iso_value_,
                      image_values_to_subdomain_indices_,
                      value_outside_),
             p::bounding_object = Mesh_3::internal::compute_bounding_box(image_),
             p::relative_error_bound = relative_error_bound_,
             p::p_rng = p_rng_,
             p::null_subdomain_index =
                     create_null_subdomain_index(null_subdomain_index_),
             p::construct_surface_patch_index =
                     create_construct_surface_patch_index(construct_surface_patch_index_));

}

template<class BGT, class Subdomain_index, class Surface_patch_index>
template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
auto
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_labeled_image_mesh_domain(const CGAL::Image_3& image_, const CGAL_NP_CLASS& np)
{
  using parameters::get_parameter;
  using parameters::get_parameter_reference;
  using parameters::choose_parameter;

  auto iso_value_ = choose_parameter(get_parameter(np, internal_np::iso_value_param), 0);
  auto value_outside_ = choose_parameter(get_parameter(np, internal_np::voxel_value), 0);
  FT relative_error_bound_ = choose_parameter(get_parameter(np, internal_np::error_bound), FT(1e-3));
  auto image_values_to_subdomain_indices_ = choose_parameter(get_parameter(np, internal_np::image_subdomain_index), Null_functor());
  CGAL::Random* p_rng_ = choose_parameter(get_parameter(np, internal_np::rng), nullptr);
  auto null_subdomain_index_ = choose_parameter(get_parameter(np, internal_np::null_subdomain_index_param), Null_functor());
  auto construct_surface_patch_index_ = choose_parameter(get_parameter(np, internal_np::surface_patch_index), Null_functor());

  using Image_ref_type = typename internal_np::Lookup_named_param_def<internal_np::weights_param_t,
                                                                      CGAL_NP_CLASS,
                                                                      CGAL::Image_3>::reference;
  CGAL::Image_3 no_weights_;
  Image_ref_type weights_ = choose_parameter(get_parameter_reference(np, internal_np::weights_param), no_weights_);
  auto features_detector_ = choose_parameter(get_parameter(np, internal_np::features_detector_param), Null_functor());

  using Default_input_features = std::vector<std::vector<typename Labeled_mesh_domain_3::Point_3>>;
  using Input_features_ref_type = typename internal_np::Lookup_named_param_def<internal_np::input_features_param_t,
                                                                               CGAL_NP_CLASS,
                                                                               Default_input_features>::reference;
  Default_input_features empty_vec;
  Input_features_ref_type input_features_
        = choose_parameter(get_parameter_reference(np, internal_np::input_features_param), empty_vec);

  CGAL_USE(iso_value_);
  namespace p = CGAL::parameters;

  const bool use_weights = weights_.is_valid();
  auto image_wrapper = use_weights
    ? create_weighted_labeled_image_wrapper(image_,
                                            weights_,
                                            image_values_to_subdomain_indices_,
                                            value_outside_)
    : create_labeled_image_wrapper(image_,
                                   image_values_to_subdomain_indices_,
                                   value_outside_);

  // warning : keep Return_type consistent with actual return type
  const bool no_features
    =  CGAL::parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::features_detector_param_t>::value
    && CGAL::parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::input_features_param_t>::value;
  using Return_type = std::conditional_t <
    no_features,
    Labeled_mesh_domain_3,
    Mesh_domain_with_polyline_features_3<Labeled_mesh_domain_3>
  >;

  Return_type domain
    (p::function = image_wrapper,
     p::bounding_object = Mesh_3::internal::compute_bounding_box(image_),
     p::relative_error_bound = relative_error_bound_,
     p::p_rng = p_rng_,
     p::null_subdomain_index =
             create_null_subdomain_index(null_subdomain_index_),
     p::construct_surface_patch_index =
             create_construct_surface_patch_index(construct_surface_patch_index_));

  // features
  Mesh_3::internal::Add_features_in_domain<!no_features>()
    (image_, weights_, domain, input_features_, features_detector_);

  return domain;
}

template<class BGT, class Subdomain_index, class Surface_patch_index>
template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_gray_image_mesh_domain(const CGAL_NP_CLASS& np)
{
  static_assert(!parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::image_3_param_t>::value, "Value for required parameter not found");
  using parameters::get_parameter;
  using parameters::get_parameter_reference;
  using parameters::choose_parameter;
  const CGAL::Image_3& image_ = get_parameter_reference(np,internal_np::image_3_param);
  auto iso_value_ = choose_parameter(get_parameter(np, internal_np::iso_value_param), 0);
  auto value_outside_ = choose_parameter(get_parameter(np, internal_np::voxel_value), 0);
  FT relative_error_bound_ = choose_parameter(get_parameter(np, internal_np::error_bound), FT(1e-3));
  auto image_values_to_subdomain_indices_ = choose_parameter(get_parameter(np, internal_np::image_subdomain_index), Null_functor());
  CGAL::Random* p_rng_ = choose_parameter(get_parameter(np, internal_np::rng), nullptr);
  auto null_subdomain_index_ = choose_parameter(get_parameter(np, internal_np::null_subdomain_index_param), Null_functor());
  auto construct_surface_patch_index_ = choose_parameter(get_parameter(np, internal_np::surface_patch_index), Null_functor());
  namespace p = CGAL::parameters;
  return Labeled_mesh_domain_3
          (p::function = create_gray_image_wrapper
                   (image_,
                    iso_value_,
                    image_values_to_subdomain_indices_,
                    value_outside_),
           p::bounding_object = Mesh_3::internal::compute_bounding_box(image_),
           p::relative_error_bound = relative_error_bound_,
           p::p_rng = p_rng_,
           p::null_subdomain_index =
                   create_null_subdomain_index(null_subdomain_index_),
           p::construct_surface_patch_index =
                   create_construct_surface_patch_index(construct_surface_patch_index_));

}

template<class BGT, class Subdomain_index, class Surface_patch_index>
template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_gray_image_mesh_domain(const CGAL::Image_3& image_,
                              const CGAL_NP_CLASS_1&  np1,
                              const CGAL_NP_CLASS_2&  np2,
                              const NP& ... nps)
{
  return create_gray_image_mesh_domain(image_, internal_np::combine_named_parameters(np1, np2, nps...));
}

template<class BGT, class Subdomain_index, class Surface_patch_index>
template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_gray_image_mesh_domain(const CGAL_NP_CLASS_1&  np1,
                              const CGAL_NP_CLASS_2&  np2,
                              const NP& ... nps)
{
  return create_gray_image_mesh_domain(internal_np::combine_named_parameters(np1, np2, nps...));
}

template<class BGT, class Subdomain_index, class Surface_patch_index>
template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
auto
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_labeled_image_mesh_domain(const CGAL_NP_CLASS& np)
{
  static_assert(!parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::image_3_param_t>::value, "Value for required parameter not found");
  using parameters::get_parameter_reference;
  const CGAL::Image_3& image_ = get_parameter_reference(np,internal_np::image_3_param);
  return create_labeled_image_mesh_domain(image_, np);
}

template<class BGT, class Subdomain_index, class Surface_patch_index>
template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
auto
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_labeled_image_mesh_domain(const CGAL::Image_3& image_,
                                 const CGAL_NP_CLASS_1&  np1,
                                 const CGAL_NP_CLASS_2&  np2,
                                 const NP& ... nps)
{
  return create_labeled_image_mesh_domain(image_, internal_np::combine_named_parameters(np1, np2, nps...));
}

template<class BGT, class Subdomain_index, class Surface_patch_index>
template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
auto
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_labeled_image_mesh_domain(const CGAL_NP_CLASS_1& np1,
                                 const CGAL_NP_CLASS_2& np2,
                                 const NP& ... nps)
{
  return create_labeled_image_mesh_domain(internal_np::combine_named_parameters(np1, np2, nps...));
}

#ifndef CGAL_NO_DEPRECATED_CODE
template<class BGT, class Subdomain_index, class Surface_patch_index>
template<typename SubdomainIndex,
         typename NullSubdomainIndex,
         typename ConstructSurfacePatchIndex>
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_gray_image_mesh_domain(const CGAL::Image_3& image_,
                              double iso_value,
                              double value_outside,
                              double relative_error_bound,
                              CGAL::Random* rng,
                              SubdomainIndex image_values_to_subdom_indices,
                              NullSubdomainIndex null_subdomain_index_,
                              ConstructSurfacePatchIndex construct_surface_patch_index_)
{
  return create_gray_image_mesh_domain(image_, parameters::iso_value(iso_value)
                                                          .image_values_to_subdomain_indices(image_values_to_subdom_indices)
                                                          .value_outside(value_outside)
                                                          .relative_error_bound(relative_error_bound)
                                                          .p_rng(rng).null_subdomain_index(null_subdomain_index_)
                                                          .construct_surface_patch_index(construct_surface_patch_index_));
}

template<class BGT, class Subdomain_index, class Surface_patch_index>
template<typename SubdomainIndex,
         typename NullSubdomainIndex,
         typename ConstructSurfacePatchIndex>
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_labeled_image_mesh_domain(const CGAL::Image_3& image_,
                                 double relative_error_bound,
                                 const CGAL::Image_3& weights_,
                                 int value_outside,
                                 CGAL::Random* rng,
                                 SubdomainIndex image_values_to_subdom_indices,
                                 NullSubdomainIndex null_subdomain_index_,
                                 ConstructSurfacePatchIndex construct_surface_patch_index_)
{
  return create_labeled_image_mesh_domain(image_, parameters::weights(weights_)
                                                             .image_values_to_subdomain_indices(image_values_to_subdom_indices)
                                                             .value_outside(value_outside)
                                                             .relative_error_bound(relative_error_bound)
                                                             .p_rng(rng)
                                                             .null_subdomain_index(null_subdomain_index_)
                                                             .construct_surface_patch_index(construct_surface_patch_index_));
}
#endif

template<class BGT, class Subdomain_index, class Surface_patch_index>
template <typename Image_word_type,
          typename FT_, typename FT2, typename Functor>
typename Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::Function_type
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_gray_image_wrapper_with_known_word_type
(const CGAL::Image_3& image ,
 const FT_& iso_value,
 const Functor& image_values_to_subdomain_indices,
 const FT2& value_outside)
{
  using Mesh_3::internal::Create_gray_image_values_to_subdomain_indices;
  typedef Create_gray_image_values_to_subdomain_indices<Functor> C_i_v_t_s_i;
  typedef typename C_i_v_t_s_i::type Image_values_to_subdomain_indices;
  Image_values_to_subdomain_indices transform_fct =
    C_i_v_t_s_i()(image_values_to_subdomain_indices, iso_value);

  typedef Mesh_3::Image_to_labeled_function_wrapper<Image_word_type,
                                                    double,
                                                    Subdomain_index,
                                                    false>           Wrapper;
  return Wrapper(image,
                 transform_fct,
                 value_outside) ;
}

template<class BGT, class Subdomain_index, class Surface_patch_index>
template <typename FT_, typename FT2, typename Functor>
typename Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::Function_type
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_gray_image_wrapper(const CGAL::Image_3& image,
                          const FT_& iso_value,
                          const Functor& image_values_to_subdomain_indices,
                          const FT2& value_outside)
{
  CGAL_IMAGE_IO_CASE(image.image(),
     return create_gray_image_wrapper_with_known_word_type<Word>
                     (image,
                      iso_value,
                      image_values_to_subdomain_indices,
                      value_outside);
                     );
  CGAL_error_msg("This place should never be reached, because it would mean "
                 "the image word type is a type that is not handled by "
                 "CGAL_ImageIO.");
  return Function_type();
}

template<class BGT, class Subdomain_index, class Surface_patch_index>
template <typename Image_word_type,
          typename FT_, typename Functor>
typename Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::Function_type
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_labeled_image_wrapper_with_known_word_type
(const CGAL::Image_3& image,
 const Functor& image_values_to_subdomain_indices,
 const FT_& value_outside)
{
  using Mesh_3::internal::Create_labeled_image_values_to_subdomain_indices;
  typedef Create_labeled_image_values_to_subdomain_indices<Functor> C_i_v_t_s_i;
  typedef typename C_i_v_t_s_i::type Image_values_to_subdomain_indices;
  Image_values_to_subdomain_indices transform_fct =
    C_i_v_t_s_i()(image_values_to_subdomain_indices);

  typedef Mesh_3::Image_to_labeled_function_wrapper<Image_word_type,
                                                    int,
                                                    Subdomain_index> Wrapper;
  return Wrapper(image,
                 transform_fct,
                 transform_fct(value_outside));
}

template<class BGT, class Subdomain_index, class Surface_patch_index>
template <typename Image_word_type,
          typename FT_, typename Functor>
typename Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::Function_type
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_weighted_labeled_image_wrapper_with_know_word_type
(const CGAL::Image_3& image,
 const CGAL::Image_3& weights,
 const Functor& image_values_to_subdomain_indices,
 const FT_& value_outside)
{
  using Mesh_3::internal::Create_labeled_image_values_to_subdomain_indices;
  typedef Create_labeled_image_values_to_subdomain_indices<Functor> C_i_v_t_s_i;
  typedef typename C_i_v_t_s_i::type Image_values_to_subdomain_indices;
  Image_values_to_subdomain_indices transform_fct =
    C_i_v_t_s_i()(image_values_to_subdomain_indices);

  typedef Mesh_3::Image_plus_weights_to_labeled_function_wrapper<
    Image_word_type,
    int, //interpolation type
    unsigned char, // Weights_type,
    Subdomain_index> Wrapper;
  return Wrapper(image,
                 weights,
                 transform_fct,
                 transform_fct(value_outside));
}

template<class BGT, class Subdomain_index, class Surface_patch_index>
template <typename FT_, typename Functor>
typename Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::Function_type
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_labeled_image_wrapper(const CGAL::Image_3& image,
                             const Functor& image_values_to_subdomain_indices,
                             const FT_& value_outside)
{
  CGAL_IMAGE_IO_CASE(image.image(),
     return create_labeled_image_wrapper_with_known_word_type<Word>
                     (image,
                      image_values_to_subdomain_indices,
                      value_outside);
                     );
  CGAL_error_msg("This place should never be reached, because it would mean "
                 "the image word type is a type that is not handled by "
                 "CGAL_ImageIO.");
  return Function_type();
}

template<class BGT, class Subdomain_index, class Surface_patch_index>
template <typename FT_, typename Functor>
typename Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::Function_type
Labeled_mesh_domain_3<BGT, Subdomain_index, Surface_patch_index>::
create_weighted_labeled_image_wrapper(const CGAL::Image_3& image,
                                      const CGAL::Image_3& weights,
                                      const Functor& image_values_to_subdomain_indices,
                                      const FT_& value_outside)
{
  CGAL_IMAGE_IO_CASE(image.image(),
    return create_weighted_labeled_image_wrapper_with_know_word_type<Word>
                      (image,
                       weights,
                       image_values_to_subdomain_indices,
                       value_outside);
                      );
  CGAL_error_msg("This place should never be reached, because it would mean "
    "the image word type is a type that is not handled by "
    "CGAL_ImageIO.");
  return Function_type();
}

} // end of CGAL

#endif // CGAL_LABELED_MESH_DOMAIN_3_IMAGE_H
