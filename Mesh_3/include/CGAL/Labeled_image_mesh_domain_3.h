// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************

#ifndef LABELED_IMAGE_MESH_DOMAIN_3_H
#define LABELED_IMAGE_MESH_DOMAIN_3_H


#include <CGAL/Mesh_3/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/Image_to_labeled_function_wrapper.h>


namespace CGAL {

/**
 * @class Labeled_image_mesh_domain_3
 *
 *
 */
template<class Image, class BGT>
class Labeled_image_mesh_domain_3
: public Mesh_3::Labeled_mesh_domain_3<
                Mesh_3::Image_to_labeled_function_wrapper<Image, BGT>,
                BGT>
{
public:
  typedef Mesh_3::Image_to_labeled_function_wrapper<Image, BGT> Image_wrapper;
  typedef Mesh_3::Labeled_mesh_domain_3<Image_wrapper, BGT> Base;

  typedef typename Base::Sphere_3 Sphere_3;
  typedef typename Base::FT FT;
  typedef BGT Geom_traits;
//  typedef typename BGT::Iso_cuboid_3 Bbox_3;
  typedef CGAL::Bbox_3 Bbox_3;

  /// Constructor
  Labeled_image_mesh_domain_3(const Image& image,
                               const FT& error_bound = FT(1e-3))
    : Base(Image_wrapper(image),
           compute_bounding_box(image),
           error_bound)                           { };

  /// Destructor
  virtual ~Labeled_image_mesh_domain_3() { };


private:
  /// Returns a sphere enclosing image \c im
  Sphere_3 compute_bounding_sphere(const Image& im) const;
  Bbox_3 compute_bounding_box(const Image& im) const
  {
    return Bbox_3(0,0,0,
                  im.xdim()*im.vx(), im.ydim()*im.vy(), im.zdim()*im.vz());
  }

private:
  // Disabled copy constructor & assignment operator
  typedef Labeled_image_mesh_domain_3<Image, BGT> Self;
  Labeled_image_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Labeled_image_mesh_domain_3





template<class Im, class BGT>
typename Labeled_image_mesh_domain_3<Im,BGT>::Sphere_3
Labeled_image_mesh_domain_3<Im,BGT>::compute_bounding_sphere(const Im& im) const
{
  typedef typename Base::Point_3 Point_3;

  // Get center of image and return a sphere centered on center with radius
  // [origin,center]
  const Point_3 center(im.xdim()*im.vx()/2.,
                       im.ydim()*im.vy()/2.,
                       im.zdim()*im.vz()/2.);

  const FT sq_radius = CGAL::squared_distance(center, Point_3(CGAL::ORIGIN));

  // We add 1 to radius to ensure that image is strictly inside sphere
  return Sphere_3(center, sq_radius + 1.);
}




}  // end namespace CGAL



#endif // LABELED_IMAGE_MESH_DOMAIN_3_H
