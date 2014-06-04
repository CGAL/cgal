/*
 * Periodic_labeled_mesh_domain_3.h
 *
 *  Created on: May 27, 2014
 *      Author: apelle
 */

#ifndef CGAL_PERIODIC_LABELED_MESH_DOMAIN_3_H
#define CGAL_PERIODIC_LABELED_MESH_DOMAIN_3_H

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>

#include <CGAL/Periodic_mesh_3/config.h>


namespace CGAL
{
template<class Function, class BGT>
class Periodic_labeled_mesh_domain_3 : public Labeled_mesh_domain_3<Function, BGT>
{
public:
  /// Base type
  typedef Labeled_mesh_domain_3<Function, BGT> Base;

  /// Public types
  typedef typename Base::FT FT;
  typedef BGT Geom_traits;
  typedef typename Base::Sphere_3 Sphere_3;
  typedef typename Base::Bbox_3 Bbox_3;
  typedef typename Base::Iso_cuboid_3 Iso_cuboid_3;

  Periodic_labeled_mesh_domain_3(const Function& f,
                         const Iso_cuboid_3& bbox,
                         const FT& error_bound = FT(1e-3))
  : Base(f, bbox, error_bound)
  {
  }

  const Iso_cuboid_3& periodic_bounding_box() const { return Base::bounding_box(); }

private:
  // Disabled copy constructor & assignment operator
  Periodic_labeled_mesh_domain_3(const Periodic_labeled_mesh_domain_3&);
  Periodic_labeled_mesh_domain_3& operator=(const Periodic_labeled_mesh_domain_3&);
};
}

#endif /* CGAL_PERIODIC_LABELED_MESH_DOMAIN_3_H */
