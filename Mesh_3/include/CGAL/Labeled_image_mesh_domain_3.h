// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************

#ifndef CGAL_LABELED_IMAGE_MESH_DOMAIN_3_H
#define CGAL_LABELED_IMAGE_MESH_DOMAIN_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Random.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/Image_to_labeled_function_wrapper.h>
#include <CGAL/Mesh_3/Null_subdomain_index.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Default.h>

namespace CGAL {

/*!
\ingroup PkgMesh3Domains

\deprecated The class template `Labeled_image_mesh_domain_3` is deprecated
since CGAL-4.13, in favor of the class template `Labeled_mesh_domain_3` and
its static function
`Labeled_mesh_domain_3::create_labeled_image_mesh_domain()`.

The class `Labeled_image_mesh_domain_3` implements a domain described by a 3D labeled image. A 3D
labeled image is a grid of voxels, where each voxel is associated with an index
(a subdomain index) characterizing the subdomain in which the voxel lies. This
class is a model of the concept `MeshDomain_3`. The domain to be discretized
is the union of voxels that have an non-default index (different from the
default constructed value of the type `Image::Type`).

This class includes a member function that provides, by interpolation, the index
of the subdomain in which any query point lies. An intersection between a segment and bounding
surfaces is detected when both segment endpoints are associated with different
values of subdomain indices. The intersection is then constructed by bisection.
The bisection stops when the query segment is shorter than a given error bound
`e`. This error bound is given by `e=d`\f$ \times\f$`bound` where `d` is the
length of the diagonal of the bounding box (in world coordinates) and
`bound` is the argument passed to the constructor of `Labeled_image_mesh_domain_3`.


\tparam Image is the type of the input image.
This parameter must be `CGAL::Image_3`.

\tparam BGT is a geometric traits class which provides
the basic operations to implement
intersection tests and intersection computations
through a bisection method. This parameter must be instantiated
with a model of the concept `BisectionGeometricTraits_3`.

\cgalModels `MeshDomain_3`

An executable that uses `Labeled_image_mesh_domain_3` must be linked with
the <I>CGAL_ImageIO</I> library.

@todo Document or comment the other parameters

\sa `BisectionGeometricTraits_3`
\sa `CGAL::make_mesh_3()`.

*/

template<class Image,
         class BGT,
         typename Image_word_type_ = unsigned char,
         typename Subdomain_index = int,
         class Null_subdomain_index = Default,
         class Wrapper_ = Default >
class
CGAL_DEPRECATED_MSG
( "The class template `CGAL::Labeled_image_mesh_domain_3` is now deprecated. "
  "Use the static member function template "
  "`Labeled_mesh_domain_3<K>::create_labeled_image_mesh_domain` instead.")
Labeled_image_mesh_domain_3
  : public Labeled_mesh_domain_3<BGT, Subdomain_index>
{
public:
  typedef Image_word_type_ Image_word_type;
  typedef typename Default::Get
    <Wrapper_,
     Mesh_3::Image_to_labeled_function_wrapper<Image_word_type,
                                               int,
                                               Subdomain_index>
     >::type Wrapper;
  typedef typename Default::Get<Null_subdomain_index,
                                CGAL::Null_subdomain_index>::type Null;

  typedef Labeled_mesh_domain_3<BGT, Subdomain_index> Base;

  typedef typename Base::Sphere_3 Sphere_3;
  typedef typename Base::FT FT;
  typedef BGT Geom_traits;
  typedef CGAL::Bbox_3 Bbox_3;
  typedef CGAL::Identity<Subdomain_index> Identity;


  /// \name Creation
  /// @{

  /*!
    Construction from an image.
    The parameter `error_bound` is relative to the size of the image.
    @todo Document or comment the other parameters
  */
  Labeled_image_mesh_domain_3(const Image& image,
                              const FT& error_bound = FT(1e-3),
                              Subdomain_index value_outside = 0,
                              Null null = Null(),
                              CGAL::Random* p_rng = nullptr)
    : Base(parameters::function = Wrapper(image, Identity(), value_outside),
           parameters::bounding_object = compute_bounding_box(image),
           parameters::relative_error_bound = error_bound,
           parameters::null_subdomain_index = null,
           parameters::p_rng = p_rng)
  {}

  /// @}

  Labeled_image_mesh_domain_3(const Image& image,
                              const FT error_bound,
                              CGAL::Random* p_rng)
    : Base(parameters::function = Wrapper(image),
           parameters::bounding_object = compute_bounding_box(image),
           parameters::relative_error_bound = error_bound,
           parameters::p_rng = p_rng)
  {}

  /// Destructor
  virtual ~Labeled_image_mesh_domain_3() {}

  using Base::bbox;

private:
  /// Returns a box enclosing image `im`
  Bbox_3 compute_bounding_box(const Image& im) const
  {
    return Bbox_3(-im.vx()+im.tx(),
                  -im.vy()+im.ty(),
                  -im.vz()+im.tz(),
                  double(im.xdim()+1)*im.vx()+im.tx(),
                  double(im.ydim()+1)*im.vy()+im.ty(),
                  double(im.zdim()+1)*im.vz()+im.tz());
  }
};  // end class Labeled_image_mesh_domain_3



}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_LABELED_IMAGE_MESH_DOMAIN_3_H
