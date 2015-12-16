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
  typedef typename Base::Surface_patch Surface_patch;
  typedef typename Base::Surface_patch_index Surface_patch_index;
  typedef typename Base::Index Index;
  typedef typename Base::Intersection Intersection;
  typedef typename Base::Subdomain_index Subdomain_index;
  typedef typename Base::Segment_3 Segment_3;
  typedef typename Base::Ray_3 Ray_3;
  typedef typename Base::Line_3 Line_3;
  typedef typename Base::Point_3 Point_3;

  Periodic_labeled_mesh_domain_3(const Function& f,
                         const Iso_cuboid_3& bbox,
                         const FT& error_bound = FT(1e-3))
  : Base(f, bbox, error_bound)
  {
  }

  const Iso_cuboid_3& periodic_bounding_box() const { return Base::bounding_box(); }

  struct Do_intersect_surface
  {
    Do_intersect_surface(const Periodic_labeled_mesh_domain_3& domain)
      : r_domain_(domain) {}

    Surface_patch operator()(const Segment_3& s) const
    {
      return this->operator()(s.source(), s.target());
    }

    Surface_patch operator()(const Ray_3& r) const
    {
      return clip_to_segment(r);
    }

    Surface_patch operator()(const Line_3& l) const
    {
      return clip_to_segment(l);
    }

  private:
    /// Returns true if points \c a & \c b do not belong to the same subdomain
    /// \c index is set to the surface index of subdomains f(a), f(b)
    Surface_patch operator()(const Point_3& a, const Point_3& b) const
    {
      // If f(a) != f(b), then [a,b] intersects some surface. Here we consider
      // [a,b] intersects surface_patch labelled <f(a),f(b)> (or <f(b),f(a)>).
      // It may be false, further rafinement will improve precision

      Iso_cuboid_3 pbb = r_domain_.periodic_bounding_box();
      FT dimension [3] = { pbb.xmax()-pbb.xmin(), pbb.ymax()-pbb.ymin(), pbb.zmax()-pbb.zmin() };
      FT a_t [3] = { a.x(), a.y(), a.z() };
      FT b_t [3] = { b.x(), b.y(), b.z() };
      a_t[0] -= pbb.xmin();
      a_t[1] -= pbb.ymin();
      a_t[2] -= pbb.zmin();
      b_t[0] -= pbb.xmin();
      b_t[1] -= pbb.ymin();
      b_t[2] -= pbb.zmin();
      int o1 [3] = { a_t[0] / dimension[0], a_t[1] / dimension[1], a_t[2] / dimension[2] };
      int o2 [3] = { b_t[0] / dimension[0], b_t[1] / dimension[1], b_t[2] / dimension[2] };

      FT a_min [3] = { a.x(), a.y(), a.z() };
      FT b_min [3] = { b.x(), b.y(), b.z() };
      for (unsigned idx = 0; idx < 3; ++idx)
      {
        FT offset = dimension[idx] * static_cast<FT>((std::min)(o1[idx], o2[idx]));
        a_min[idx] -= offset;
        b_min[idx] -= offset;
      }

      Subdomain_index value_a = r_domain_.labeling_function()(Point_3(a_min[0], a_min[1], a_min[2]));
      Subdomain_index value_b = r_domain_.labeling_function()(Point_3(b_min[0], b_min[1], b_min[2]));
      if ( value_a != value_b )
        return Surface_patch(r_domain_.make_surface_index(value_a, value_b));

      FT a_max [3] = { a.x(), a.y(), a.z() };
      FT b_max [3] = { b.x(), b.y(), b.z() };
      for (unsigned idx = 0; idx < 3; ++idx)
      {
        FT offset = dimension[idx] * static_cast<FT>((std::max)(o1[idx], o2[idx]));
        a_max[idx] -= offset;
        b_max[idx] -= offset;
      }

      value_a = r_domain_.labeling_function()(Point_3(a_max[0], a_max[1], a_max[2]));
      value_b = r_domain_.labeling_function()(Point_3(b_max[0], b_max[1], b_max[2]));
      if ( value_a != value_b )
        return Surface_patch(r_domain_.make_surface_index(value_a, value_b));

      return Surface_patch();
    }

    /**
     * Clips \c query to a segment \c s, and call operator()(s)
     */
    template<typename Query>
    Surface_patch clip_to_segment(const Query& query) const
    {
      typename cpp11::result_of<typename BGT::Intersect_3(Query, Iso_cuboid_3)>::type
        clipped = CGAL::intersection(query, r_domain_.bounding_box());

      if(clipped)
#if CGAL_INTERSECTION_VERSION > 1
        if(const Segment_3* s = boost::get<Segment_3>(&*clipped))
          return this->operator()(*s);
#else
        if(const Segment_3* s = object_cast<Segment_3>(&clipped))
          return this->operator()(*s);
#endif

      return Surface_patch();
    }

  private:
    const Periodic_labeled_mesh_domain_3& r_domain_;
  };

  /// Returns Do_intersect_surface object
  Do_intersect_surface do_intersect_surface_object() const
  {
    return Do_intersect_surface(*this);
  }

  struct Construct_intersection
  {
    Construct_intersection(const Periodic_labeled_mesh_domain_3& domain)
      : r_domain_(domain) {}

    Intersection operator()(const Segment_3& s) const
    {
#ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
      CGAL_precondition(r_domain_.do_intersect_surface_object()(s));
#endif // NOT CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
      return this->operator()(s.source(),s.target());
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
    /**
     * Returns a point in the intersection of [a,b] with the surface
     * \c a must be the source point, and \c b the out point. It's important
     * because it drives bisection cuts.
     * Indeed, the returned point is the first intersection from \c [a,b]
     * with a subdomain surface.
     */
    Intersection operator()(const Point_3& a, const Point_3& b) const
    {
      // Functors
      typename BGT::Compute_squared_distance_3 squared_distance =
                                      BGT().compute_squared_distance_3_object();
      typename BGT::Construct_midpoint_3 midpoint =
                                      BGT().construct_midpoint_3_object();

      Iso_cuboid_3 pbb = r_domain_.periodic_bounding_box();
      FT dimension [3] = { pbb.xmax()-pbb.xmin(), pbb.ymax()-pbb.ymin(), pbb.zmax()-pbb.zmin() };
      FT a_t [3] = { a.x(), a.y(), a.z() };
      FT b_t [3] = { b.x(), b.y(), b.z() };
      a_t[0] -= pbb.xmin();
      a_t[1] -= pbb.ymin();
      a_t[2] -= pbb.zmin();
      b_t[0] -= pbb.xmin();
      b_t[1] -= pbb.ymin();
      b_t[2] -= pbb.zmin();

      int o1 [3] = { a_t[0] / dimension[0], a_t[1] / dimension[1], a_t[2] / dimension[2] };
      int o2 [3] = { b_t[0] / dimension[0], b_t[1] / dimension[1], b_t[2] / dimension[2] };

      FT a_min [3] = { a.x(), a.y(), a.z() };
      FT b_min [3] = { b.x(), b.y(), b.z() };
      for (unsigned idx = 0; idx < 3; ++idx)
      {
        FT offset = dimension[idx] * static_cast<FT>((std::min)(o1[idx], o2[idx]));
        a_min[idx] -= offset;
        b_min[idx] -= offset;
      }

      // Non const points
      Point_3 p1(a_min[0], a_min[1], a_min[2]);
      Point_3 p2(b_min[0], b_min[1], b_min[2]);
      Point_3 mid = midpoint(p1, p2);

      // Cannot be const: those values are modified below.
      Subdomain_index value_at_p1 = r_domain_.labeling_function()(p1);
      Subdomain_index value_at_p2 = r_domain_.labeling_function()(p2);
      Subdomain_index value_at_mid = r_domain_.labeling_function()(mid,true);

      // If both extremities are in the same subdomain,
      // there is no intersection.
      // This should not happen...
      if( value_at_p1 == value_at_p2 )
      {
        FT a_max [3] = { a.x(), a.y(), a.z() };
        FT b_max [3] = { b.x(), b.y(), b.z() };
        for (unsigned idx = 0; idx < 3; ++idx)
        {
          FT offset = dimension[idx] * static_cast<FT>((std::max)(o1[idx], o2[idx]));
          a_max[idx] -= offset;
          b_max[idx] -= offset;
        }
        p1 = Point_3(a_max[0], a_max[1], a_max[2]);
        p2 = Point_3(b_max[0], b_max[1], b_max[2]);
        mid = midpoint(p1, p2);

        value_at_p1 = r_domain_.labeling_function()(p1);
        value_at_p2 = r_domain_.labeling_function()(p2);
        value_at_mid = r_domain_.labeling_function()(mid,true);

        if( value_at_p1 == value_at_p2 )
          return Intersection();
      }

      // Construct the surface patch index and index from the values at 'a'
      // and 'b'. Even if the bissection find out a different pair of
      // values, the reported index will be constructed from the initial
      // values.
      const Surface_patch_index sp_index =
        r_domain_.make_surface_index(value_at_p1, value_at_p2);
      const Index index = r_domain_.index_from_surface_patch_index(sp_index);

      // Else lets find a point (by bisection)
      // Bisection ends when the point is near than error bound from surface
      while(true)
      {
        // If the two points are enough close, then we return midpoint
        if ( squared_distance(p1, p2) < r_domain_.squared_error_bound_value() )
        {
          CGAL_assertion(value_at_p1 != value_at_p2);
          return Intersection(mid, index, 2);
        }

        // Else we must go on
        // Here we consider that p1(a) is the source point. Thus, we keep p1 and
        // change p2 if f(p1)!=f(p2).
        // That allows us to find the first intersection from a of [a,b] with
        // a surface.
        if ( value_at_p1 != value_at_mid )
        {
          p2 = mid;
          value_at_p2 = value_at_mid;
        }
        else
        {
          p1 = mid;
          value_at_p1 = value_at_mid;
        }

        mid = midpoint(p1, p2);
        value_at_mid = r_domain_.labeling_function()(mid,true);
      }
    }

    /// Clips \c query to a segment \c s, and call operator()(s)
    template<typename Query>
    Intersection clip_to_segment(const Query& query) const
    {
      typename cpp11::result_of<typename BGT::Intersect_3(Query, Iso_cuboid_3)>::type
        clipped = CGAL::intersection(query, r_domain_.bounding_box());

      if(clipped)
#if CGAL_INTERSECTION_VERSION > 1
        if(const Segment_3* s = boost::get<Segment_3>(&*clipped))
          return this->operator()(*s);
#else
        if(const Segment_3* s = object_cast<Segment_3>(&clipped))
          return this->operator()(*s);
#endif

      return Intersection();
    }

  private:
    const Periodic_labeled_mesh_domain_3& r_domain_;
  };

  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection(*this);
  }

private:
  // Disabled copy constructor & assignment operator
  Periodic_labeled_mesh_domain_3(const Periodic_labeled_mesh_domain_3&);
  Periodic_labeled_mesh_domain_3& operator=(const Periodic_labeled_mesh_domain_3&);
};
}

#endif /* CGAL_PERIODIC_LABELED_MESH_DOMAIN_3_H */
