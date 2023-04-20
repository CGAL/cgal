// Copyright (c) 2021 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Mael Rouxel-Labbé

#ifndef CGAL_AABB_TREE_INTERNAL_TRIANGLE_DATUM_COVERING_H
#define CGAL_AABB_TREE_INTERNAL_TRIANGLE_DATUM_COVERING_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/AABB_tree/internal/AABB_traversal_traits.h>
#include <CGAL/AABB_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/array.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Container_helper.h>
#include <CGAL/property_map.h>

#include <queue>
#include <unordered_set>
#include <utility>
#include <vector>

namespace CGAL {
namespace AABB_trees {
namespace internal {

// std::vector to property map
// Not using Pointer_property_map because the underlying range is empty initially and can change
template <typename T>
struct Vector_property_map
{
  using Range = std::vector<T>;

  using key_type = std::size_t;
  using value_type = T;
  using reference = value_type&;
  using category = boost::read_write_property_map_tag;

  Vector_property_map() : m_range_ptr(std::make_shared<Range>()) { }

  inline friend void put(const Vector_property_map& map, const key_type& k, const value_type& v)
  {
    CGAL_precondition(map.m_range_ptr != nullptr);

    if(k >= map.m_range_ptr->size())
      map.m_range_ptr->resize(k+1);

    map.m_range_ptr->operator[](k) = v;
  }

  inline friend reference get(const Vector_property_map& map, const key_type& k)
  {
    CGAL_precondition(map.m_range_ptr != nullptr);
    return map.m_range_ptr->operator[](k);
  }

  Range& range() { return *m_range_ptr; }
  const Range& range() const { return *m_range_ptr; }

private:
  std::shared_ptr<Range> m_range_ptr;
};

// Same as the standard traversal traits, but for multiple primitives per datum,
// such that the final operation on the datum is only performed once.
template <typename AABBTraits, typename BaseTraversalTraits>
struct Covered_traversal_traits
  : BaseTraversalTraits
{
  using Base = BaseTraversalTraits;
  using Primitive = typename AABBTraits::Primitive;

  std::unordered_set<std::size_t> visited_data;

public:
  template <typename ... Args>
  Covered_traversal_traits(Args&&... args) : Base(std::forward<Args>(args)...) { }

  template <typename Query>
  void intersection(const Query& query, const Primitive& primitive)
  {
    // check a datum only once
    auto is_insert_successful = visited_data.insert(primitive.id().second/*unique input face ID*/);
    if(!is_insert_successful.second)
      return;

    return Base::intersection(query, primitive);
  }
};

template <typename AABBTraits>
struct Covered_tree_traversal_traits
{
  using Base_projection_traits = CGAL::internal::AABB_tree::Projection_traits<AABBTraits>;
  using Projection_traits = Covered_traversal_traits<AABBTraits, Base_projection_traits>;

  template <typename Query>
  using Do_intersect_traits_base = CGAL::internal::AABB_tree::Do_intersect_traits<AABBTraits, Query>;
  template <typename Query>
  using Do_intersect_traits = Covered_traversal_traits<AABBTraits, Do_intersect_traits_base<Query> >;

  template <typename Query>
  using First_intersection_traits_base = CGAL::internal::AABB_tree::First_intersection_traits<AABBTraits, Query>;
  template <typename Query>
  using First_intersection_traits = Covered_traversal_traits<AABBTraits, First_intersection_traits_base<Query> >;

  template <typename Query>
  using First_primitive_traits_base = CGAL::internal::AABB_tree::First_primitive_traits<AABBTraits, Query>;
  template <typename Query>
  using First_primitive_traits = Covered_traversal_traits<AABBTraits, First_primitive_traits_base<Query> >;

  template <typename Query, typename CountingIterator>
  using Listing_primitive_traits_base = CGAL::internal::AABB_tree::Listing_primitive_traits<AABBTraits, Query, CountingIterator>;
  template <typename Query, typename CountingIterator>
  using Listing_primitive_traits = Covered_traversal_traits<AABBTraits, Listing_primitive_traits_base<Query, CountingIterator> >;

  template <typename Query, typename CountingIterator>
  using Listing_intersection_traits_base = CGAL::internal::AABB_tree::Listing_intersection_traits<AABBTraits, Query, CountingIterator>;
  template <typename Query, typename CountingIterator>
  using Listing_intersection_traits = Covered_traversal_traits<AABBTraits, Listing_intersection_traits_base<Query, CountingIterator> >;
};

// Dissociated from the class `AABB_covered_triangle_tree` for clarity
template <typename Kernel, typename Point>
struct AABB_covered_triangle_tree_traits
{
  using Triangle_3 = typename Kernel::Triangle_3;

  // Below is a lot of trouble to cover a single datum with multiple primitives using smaller bboxes
  using ID = std::pair<std::size_t /*primitive ID*/, std::size_t /*input face ID*/>;
  using IDPM = CGAL::First_of_pair_property_map<ID>;

  // Primitive ID --> box vector pos --> Bounding Box
  using BPMB = internal::Vector_property_map<CGAL::Bbox_3>;
  using BPM = CGAL::Compose_property_map<IDPM, BPMB>;

  // Primitive ID --> point vector pos --> Reference Point
  using RPPMB = internal::Vector_property_map<Point>;
  using RPPM = CGAL::Compose_property_map<IDPM, RPPMB>;

  // Primitive ID --> Datum pos vector pos --> Datum pos --> Datum
  // The vector of data has size nf, but the vector of datum pos has size tree.size()
  using DPPMB = internal::Vector_property_map<std::size_t>; // pos --> Datum pos
  using DPPM = CGAL::Compose_property_map<IDPM, DPPMB>; // PID --> Datum pos
  using DPMB = internal::Vector_property_map<Triangle_3>; // Datum pos --> Datum
  using DPM = CGAL::Compose_property_map<DPPM, DPMB>; // PID --> Datum

  using Primitive = CGAL::AABB_primitive<ID, DPM, RPPM,
                                         CGAL::Tag_true /*external pmaps*/,
                                         CGAL::Tag_false /*no caching*/>;

  using AABB_geom_traits = Kernel;
  using AABB_traits = CGAL::AABB_traits<AABB_geom_traits, Primitive, BPM>;
  using AABB_tree = CGAL::AABB_tree<AABB_traits>;
};

template <typename Kernel,
          typename Point = typename Kernel::Point_3>
struct AABB_covered_triangle_tree
  : public AABB_covered_triangle_tree_traits<Kernel, Point>::AABB_tree
{
  using FT = typename Kernel::FT;
  using Triangle_3 = typename Kernel::Triangle_3;

  using ACTTT = AABB_covered_triangle_tree_traits<Kernel, Point>;

  using BPM = typename ACTTT::BPM;
  using RPPM = typename ACTTT::RPPM;
  using DPPMB = typename ACTTT::DPPMB;
  using DPPM = typename ACTTT::DPPM;
  using DPMB = typename ACTTT::DPMB;
  using DPM = typename ACTTT::DPM;

  using ID = typename ACTTT::ID;
  using Primitive = typename ACTTT::Primitive;
  using AABB_traits = typename ACTTT::AABB_traits;
  using AABB_tree = typename ACTTT::AABB_tree;
  using Base = AABB_tree;

protected:
  double m_sq_length;

  DPPMB m_dppmb; // std::size_t (id) --> datum pos

  BPM m_bpm; // std::size_t (id) --> bounding box
  RPPM m_rppm; // std::size_t (id) --> reference point
  DPMB m_dpmb; // std::size_t (datum pos) --> triangle datum

  DPM m_dpm; // std::size_t (id) --> triangle (datum)

  std::size_t fid = 0;

public:
  AABB_covered_triangle_tree(const double max_length,
                             const AABB_traits& traits = AABB_traits())
    : Base(traits),
      m_sq_length(square(max_length)),
      m_dppmb(), m_bpm(), m_rppm(), m_dpmb(),
      m_dpm(DPPM(Default(), m_dppmb/*first binder's value_map*/)/*second binder's key map*/, m_dpmb)
  {
    initialize_tree_property_maps();
  }

private:
  void initialize_tree_property_maps() const
  {
    // Can't be set in the default constructed traits that are passed to the base
    // since m_bpm is a member of the derived class.
    //
    // 'const_cast' because CGAL::AABB_tree only gives a const& to its traits...
    const_cast<AABB_traits&>(this->traits()).bbm = m_bpm;
    const_cast<AABB_traits&>(this->traits()).set_shared_data(m_dpm, m_rppm);
  }

  template <typename P> // Kernel is Simple_Cartesian<Interval>
  CGAL::Bbox_3 compute_bbox(const P& ap0, const P& ap1, const P& ap2)
  {
    double xmin = (CGAL::min)(ap0.x().inf(), (CGAL::min)(ap1.x().inf(), ap2.x().inf()));
    double ymin = (CGAL::min)(ap0.y().inf(), (CGAL::min)(ap1.y().inf(), ap2.y().inf()));
    double zmin = (CGAL::min)(ap0.z().inf(), (CGAL::min)(ap1.z().inf(), ap2.z().inf()));

    double xmax = (CGAL::max)(ap0.x().sup(), (CGAL::max)(ap1.x().sup(), ap2.x().sup()));
    double ymax = (CGAL::max)(ap0.y().sup(), (CGAL::max)(ap1.y().sup(), ap2.y().sup()));
    double zmax = (CGAL::max)(ap0.z().sup(), (CGAL::max)(ap1.z().sup(), ap2.z().sup()));

    return CGAL::Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
  }

  template <typename P> // Kernel is Simple_Cartesian<Interval>
  const P& compute_reference_point(const P&, const P&, const P& ap2)
  {
    return ap2; // ap2 is the midpoint when splitting
  }

public:
  void reserve(std::size_t nf)
  {
    CGAL::internal::reserve(m_dpmb.range(), m_dpmb.range().size() + nf);

    // Due to splitting, these might need more than 'nf'
    CGAL::internal::reserve(m_dppmb.range(), m_dppmb.range().size() + nf);
    CGAL::internal::reserve(m_rppm.value_map.range(), m_rppm.value_map.range().size() + nf);
    CGAL::internal::reserve(m_bpm.value_map.range(), m_bpm.value_map.range().size() + nf);
  }

  void split_and_insert(const Triangle_3& tr)
  {
    // Convert to intervals to ensure that the bounding box fully covers the subdividing triangle
    using IT = CGAL::Interval_nt<true>;
    using NT = typename IT::value_type;
    using AK = CGAL::Simple_cartesian<IT>;
    using K2AK = CGAL::Cartesian_converter<Kernel, AK>;
    using AFT = AK::FT;
    using APoint_3 = AK::Point_3;

    using APL = std::pair<APoint_3, NT>; // point and upper bound of the length of the opposite edge
    using AT = std::array<APL, 3>;

    const std::size_t data_size = m_dpmb.range().size();
    put(m_dpmb, data_size, tr);

    auto vertex = Kernel().construct_vertex_3_object();
    const Point& p0 = vertex(tr, 0);
    const Point& p1 = vertex(tr, 1);
    const Point& p2 = vertex(tr, 2);

    if(m_sq_length == FT(0)) // no splits
    {
      const std::size_t pid = this->size();
      ID id = std::make_pair(pid, fid++);

      put(m_dppmb, pid, data_size);
      put(m_rppm, id, p1); // the ref point that `One_point_from_face_descriptor_map` would give
      put(m_bpm, id, Kernel().construct_bbox_3_object()(tr));

//      std::cout << "Primitive[" << id.first << " " << id.second << "]; "
//                    << "Bbox: [" << get(m_bpm, id) << "] "
//                    << "Point: (" << get(m_rppm, id) << ") "
//                    << "Datum: [" << get(m_bpm, id) << "]" << std::endl;

      Primitive p(id/*, m_dpm, m_rppm*/); // pmaps are external, shared data
      this->insert(p);
      return;
    }

    K2AK k2ak;
    std::queue<AT> to_treat;

    const APoint_3 ap0 = k2ak(p0);
    const APoint_3 ap1 = k2ak(p1);
    const APoint_3 ap2 = k2ak(p2);
    const AFT sq_l0 = CGAL::squared_distance(ap1, ap2);
    const AFT sq_l1 = CGAL::squared_distance(ap2, ap0);
    const AFT sq_l2 = CGAL::squared_distance(ap0, ap1);

    to_treat.push(CGAL::make_array(std::make_pair(ap0, sq_l0.sup()),
                                   std::make_pair(ap1, sq_l1.sup()),
                                   std::make_pair(ap2, sq_l2.sup())));

    while(!to_treat.empty())
    {
      const AT t = std::move(to_treat.front());
      to_treat.pop();

      const APL& apl0 = t[0];
      const APL& apl1 = t[1];
      const APL& apl2 = t[2];

      int i = (apl0.second >= apl1.second)
                ? (apl0.second >= apl2.second) ? 0 : 2
                : (apl1.second >= apl2.second) ? 1 : 2;

      const NT max_sql = t[i].second;

      // If the face is too big, do a split (two small bboxes rather than a big one)
      if(max_sql > m_sq_length)
      {
        // Could be factorized, but this is simpler to read
        if(i == 0)
        {
          // 0 1 2 into 0 1 A and 0 A 2
          const APoint_3 amp = CGAL::midpoint(apl1.first, apl2.first);
          const NT sq_half_length = apl0.second / NT(4);
          const NT sq_diag_length = CGAL::squared_distance(amp, apl0.first).sup();

          to_treat.push(CGAL::make_array(std::make_pair(apl0.first, sq_half_length),
                                         std::make_pair(apl1.first, sq_diag_length),
                                         std::make_pair(amp, apl2.second)));
          to_treat.push(CGAL::make_array(std::make_pair(apl2.first, sq_diag_length),
                                         std::make_pair(apl0.first, sq_half_length),
                                         std::make_pair(amp, apl1.second)));
        }
        else if(i == 1)
        {
          // 0 1 2 into 0 1 A and 1 2 A
          const APoint_3 amp = CGAL::midpoint(apl2.first, apl0.first);
          const NT sq_half_length = apl1.second / NT(4);
          const NT sq_diag_length = CGAL::squared_distance(amp, apl1.first).sup();

          to_treat.push(CGAL::make_array(std::make_pair(apl0.first, sq_diag_length),
                                         std::make_pair(apl1.first, sq_half_length),
                                         std::make_pair(amp, apl2.second)));
          to_treat.push(CGAL::make_array(std::make_pair(apl1.first, sq_half_length),
                                         std::make_pair(apl2.first, sq_diag_length),
                                         std::make_pair(amp, apl0.second)));
        }
        else // i == 2
        {
          // 0 1 2 into 0 A 2 and 2 A 1
          const APoint_3 amp = CGAL::midpoint(apl0.first, apl1.first);
          const NT sq_half_length = apl2.second / NT(4);
          const NT sq_diag_length = CGAL::squared_distance(amp, apl2.first).sup();

          to_treat.push(CGAL::make_array(std::make_pair(apl2.first, sq_half_length),
                                         std::make_pair(apl0.first, sq_diag_length),
                                         std::make_pair(amp, apl1.second)));
          to_treat.push(CGAL::make_array(std::make_pair(apl1.first, sq_diag_length),
                                         std::make_pair(apl2.first, sq_half_length),
                                         std::make_pair(amp, apl0.second)));
        }
      }
      else // all edges have length below the threshold, create a primitive
      {
        const std::size_t pid = this->size();
        ID id = std::make_pair(pid, fid);

        put(m_dppmb, pid, data_size);

        // this is basically to_double() of an APoint_3, but FT is not necessarily 'double'
        const APoint_3& apt = compute_reference_point(apl0.first, apl1.first, apl2.first);
        put(m_rppm, id, Point((FT(apt.x().sup()) + FT(apt.x().inf())) / FT(2),
                              (FT(apt.y().sup()) + FT(apt.y().inf())) / FT(2),
                              (FT(apt.z().sup()) + FT(apt.z().inf())) / FT(2)));

        put(m_bpm, id, compute_bbox(apl0.first, apl1.first, apl2.first));

//        std::cout << "Primitive[" << std::get<0>(id) << " " << std::get<1>(id) << " " << std::get<2>(id) << "]; "
//                  << "Bbox: [" << get(m_bpm, id) << "] "
//                  << "Point: (" << get(m_rppm, id) << ") "
//                  << "Datum: [" << get(m_dpm, id) << "]" << std::endl;

        Primitive p(id/*, m_dpm, m_rppm*/); // pmaps are external, shared data
        this->insert(p);
      }
    }

    ++fid;
  }
};

} // namespace internal
} // namespace AABB_trees
} // namespace CGAL

#endif // CGAL_AABB_TREE_INTERNAL_TRIANGLE_DATUM_COVERING_H
