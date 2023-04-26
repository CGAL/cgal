// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_SPLITTING_HELPER_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_SPLITTING_HELPER_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_AABB_geom_traits.h>

#include <CGAL/AABB_tree/internal/AABB_traversal_traits.h>
#include <CGAL/AABB_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/array.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Container_helper.h>
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian_converter.h>

#include <queue>
#include <unordered_set>
#include <utility>
#include <vector>

namespace CGAL {
namespace Alpha_wraps_3 {
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

// Same as the standard traversal traits, but for multiple primitives per datum: the final operation
// done on the datum is only performed once.
template <typename AABBTraits>
struct Splitter_traversal_traits
{
  // -----------------------------------------------------------------------------------------------
  class Projection_traits
    : public CGAL::internal::AABB_tree::Projection_traits<AABBTraits>
  {
    using Base = CGAL::internal::AABB_tree::Projection_traits<AABBTraits>;
    using Point = typename AABBTraits::Point_3;
    using Primitive = typename AABBTraits::Primitive;

    std::unordered_set<std::size_t> visited_data;

public:
    Projection_traits(const Point& hint,
                      const typename Primitive::Id& hint_primitive,
                      const AABBTraits& traits)
      : Base(hint, hint_primitive, traits)
    { }

    void intersection(const Point& query, const Primitive& primitive)
    {
      // check a datum only once
      auto is_insert_successful = visited_data.insert(primitive.id().second/*unique input face ID*/);
      if(!is_insert_successful.second)
        return;

      return Base::intersection(query, primitive);
    }
  };

  // -----------------------------------------------------------------------------------------------
  template <typename Query>
  class Do_intersect_traits
    : public CGAL::internal::AABB_tree::Do_intersect_traits<AABBTraits, Query>
  {
    using Base = CGAL::internal::AABB_tree::Do_intersect_traits<AABBTraits, Query>;
    using Primitive = typename AABBTraits::Primitive;

    std::unordered_set<std::size_t> visited_data;

public:
    Do_intersect_traits(const AABBTraits& traits) : Base(traits) { }

    void intersection(const Query& query, const Primitive& primitive)
    {
      // check a datum only once
      auto is_insert_successful = visited_data.insert(primitive.id().second/*unique input face ID*/);
      if(!is_insert_successful.second)
        return;

      return Base::intersection(query, primitive);
    }
  };

  // -----------------------------------------------------------------------------------------------
  template <typename Query>
  class First_intersection_traits
    : public CGAL::internal::AABB_tree::First_intersection_traits<AABBTraits, Query>
  {
    using Base = CGAL::internal::AABB_tree::First_intersection_traits<AABBTraits, Query>;
    using Primitive = typename AABBTraits::Primitive;

    std::unordered_set<std::size_t> visited_data;

public:
    First_intersection_traits(const AABBTraits& traits) : Base(traits) { }

    void intersection(const Query& query, const Primitive& primitive)
    {
      // check a datum only once
      auto is_insert_successful = visited_data.insert(primitive.id().second/*unique input face ID*/);
      if(!is_insert_successful.second)
        return;

      return Base::intersection(query, primitive);
    }
  };
};

// Dissociated from the class `AABB_tree_oracle_splitter` for clarity (the AABB_tree is a template
// parameter of the other base of the oracle too)
template <typename Point, typename GT>
struct AABB_tree_splitter_traits
{
  using Triangle_3 = typename GT::Triangle_3;

  // Below is a lot of trouble to cover a single datum with multiple primitives using smaller bboxes

  // The input face ID serves when traversing the tree, to avoid doing the same intersection()
  // on the same datum seen from different primitives.
  //
  // Technically, FPM could type-erase the mesh and the VPM, as it currently forces all independent
  // inputs to have the same types. This is not such much of an issue for the mesh type,
  // but it can be annoying for the VPM type.
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

  using AABB_traits = CGAL::AABB_traits<GT, Primitive, BPM>;
  using AABB_tree = CGAL::AABB_tree<AABB_traits>;
};

template <bool subdivide, typename Point, typename GT>
struct AABB_tree_oracle_splitter
{
  using FT = typename GT::FT;
  using Triangle_3 = typename GT::Triangle_3;

  using ATST = AABB_tree_splitter_traits<Point, GT>;

  using BPM = typename ATST::BPM;
  using RPPM = typename ATST::RPPM;
  using DPPMB = typename ATST::DPPMB;
  using DPPM = typename ATST::DPPM;
  using DPMB = typename ATST::DPMB;
  using DPM = typename ATST::DPM;

  using ID = typename ATST::ID;
  using Primitive = typename ATST::Primitive;
  using AABB_traits = typename ATST::AABB_traits;
  using AABB_tree = typename ATST::AABB_tree;

protected:
  double m_sq_alpha;

  // one per face
  DPPMB m_dppmb; // std::size_t (id) --> datum pos

  // possibly more than one per face
  BPM m_bpm; // std::size_t (id) --> bounding box
  RPPM m_rppm; // std::size_t (id) --> reference point
  DPMB m_dpmb; // std::size_t (datum pos) --> triangle datum

  DPM m_dpm; // std::size_t (id) --> triangle (datum)

  std::size_t fid = 0;

public:
  AABB_tree_oracle_splitter(const double alpha)
    :
      m_sq_alpha(square(alpha)),
      m_dppmb(), m_bpm(), m_rppm(), m_dpmb(),
      m_dpm(DPPM(Default(), m_dppmb/*first binder's value_map*/)/*second binder's key map*/, m_dpmb)
  { }

protected:
  void initialize_tree_property_maps(const AABB_tree& tree) const
  {
    // Can't be set in the default constructed traits that are passed to the base
    // since m_bpm is a member of the derived class.
    //
    // 'const_cast' because CGAL::AABB_tree only gives a const& to its traits.
    const_cast<AABB_traits&>(tree.traits()).bbm = m_bpm;

    const_cast<AABB_traits&>(tree.traits()).set_shared_data(m_dpm, m_rppm);
  }

  void reserve(std::size_t nf)
  {
    CGAL::internal::reserve(m_dpmb.range(), m_dpmb.range().size() + nf);

    // Due to splitting, these might need more than 'nf'
    CGAL::internal::reserve(m_dppmb.range(), m_dppmb.range().size() + nf);
    CGAL::internal::reserve(m_rppm.value_map.range(), m_rppm.value_map.range().size() + nf);
    CGAL::internal::reserve(m_bpm.value_map.range(), m_bpm.value_map.range().size() + nf);
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

  void split_and_insert_datum(const Triangle_3& tr,
                              AABB_tree& tree,
                              const GT& gt)
  {
    // Convert to intervals to ensure that the bounding box fully covers the subdividing triangle
    using AK = CGAL::Simple_cartesian<CGAL::Interval_nt<true> >;
    using K2AK = CGAL::Cartesian_converter<GT, AK>;
    using AK2K = CGAL::Cartesian_converter<AK, GT>;
    using AFT = AK::FT;
    using APoint_3 = AK::Point_3;

    using APL = std::pair<APoint_3, AFT>; // point and length of the opposite edge
    using AT = std::array<APL, 3>;

    const std::size_t data_size = m_dpmb.range().size();
    put(m_dpmb, data_size, tr);

    auto vertex = gt.construct_vertex_3_object();

    const Point& p0 = vertex(tr, 0);
    const Point& p1 = vertex(tr, 1);
    const Point& p2 = vertex(tr, 2);

    if(!subdivide || m_sq_alpha == 0.) // no splits
    {
      const std::size_t pid = tree.size();
      ID id = std::make_pair(pid, fid++);

      put(m_dppmb, pid, data_size);
      put(m_rppm, id, p1); // the ref point that `One_point_from_face_descriptor_map` would give
      put(m_bpm, id, gt.construct_bbox_3_object()(tr));

//      std::cout << "Primitive[" << id.first << " " << id.second << "]; "
//                    << "Bbox: [" << get(m_bpm, id) << "] "
//                    << "Point: (" << get(m_rppm, id) << ") "
//                    << "Datum: [" << get(m_bpm, id) << "]" << std::endl;

      Primitive p(id/*, m_dpm, m_rppm*/); // pmaps are external, shared data
      tree.insert(p);
      return;
    }

    K2AK k2ak;
    AK2K ak2k;
    std::queue<AT> to_treat;

    const APoint_3 ap0 = k2ak(p0);
    const APoint_3 ap1 = k2ak(p1);
    const APoint_3 ap2 = k2ak(p2);
    const AFT sq_l0 = CGAL::squared_distance(ap1, ap2);
    const AFT sq_l1 = CGAL::squared_distance(ap2, ap0);
    const AFT sq_l2 = CGAL::squared_distance(ap0, ap1);

    to_treat.push(CGAL::make_array(std::make_pair(ap0, sq_l0),
                                   std::make_pair(ap1, sq_l1),
                                   std::make_pair(ap2, sq_l2)));

    while(!to_treat.empty())
    {
      const AT t = std::move(to_treat.front());
      to_treat.pop();

      const APL& apl0 = t[0];
      const APL& apl1 = t[1];
      const APL& apl2 = t[2];

      int i = (apl0.second.sup() >= apl1.second.sup())
                ? (apl0.second.sup() >= apl2.second.sup()) ? 0 : 2
                : (apl1.second.sup() >= apl2.second.sup()) ? 1 : 2;

      const FT max_sql = t[i].second.sup();

      // The '3 * alpha' is some empirically-determined value
      const FT sq_bound = 9. * m_sq_alpha;

      // If the face is too big, do a fake split (two small bboxes rather than a big one)
      if(max_sql > sq_bound)
      {
        // Could be factorized, but this is simpler to read
        if(i == 0)
        {
          // 0 1 2 into 0 1 A and 0 A 2
          const APoint_3 amp = CGAL::midpoint(apl1.first, apl2.first);
          to_treat.push(CGAL::make_array(std::make_pair(apl0.first, CGAL::squared_distance(apl1.first, amp)),
                                         std::make_pair(apl1.first, CGAL::squared_distance(amp, apl0.first)),
                                         std::make_pair(amp, apl2.second)));
          to_treat.push(CGAL::make_array(std::make_pair(apl2.first, CGAL::squared_distance(apl0.first, amp)),
                                         std::make_pair(apl0.first, CGAL::squared_distance(amp, apl2.first)),
                                         std::make_pair(amp, apl1.second)));
        }
        else if(i == 1)
        {
          // 0 1 2 into 0 1 A and 1 2 A
          const APoint_3 amp = CGAL::midpoint(apl2.first, apl0.first);
          to_treat.push(CGAL::make_array(std::make_pair(apl0.first, CGAL::squared_distance(apl1.first, amp)),
                                         std::make_pair(apl1.first, CGAL::squared_distance(amp, apl0.first)),
                                         std::make_pair(amp, apl2.second)));
          to_treat.push(CGAL::make_array(std::make_pair(apl1.first, CGAL::squared_distance(apl2.first, amp)),
                                         std::make_pair(apl2.first, CGAL::squared_distance(amp, apl1.first)),
                                         std::make_pair(amp, apl0.second)));
        }
        else // i == 2
        {
          // 0 1 2 into 0 A 2 and 2 A 1
          const APoint_3 amp = CGAL::midpoint(apl0.first, apl1.first);
          to_treat.push(CGAL::make_array(std::make_pair(apl2.first, CGAL::squared_distance(apl0.first, amp)),
                                         std::make_pair(apl0.first, CGAL::squared_distance(amp, apl2.first)),
                                         std::make_pair(amp, apl1.second)));
          to_treat.push(CGAL::make_array(std::make_pair(apl1.first, CGAL::squared_distance(apl2.first, amp)),
                                         std::make_pair(apl2.first, CGAL::squared_distance(amp, apl1.first)),
                                         std::make_pair(amp, apl0.second)));
        }
      }
      else // all edges have length below the threshold, create a primitive
      {
        const std::size_t pid = tree.size();
        ID id = std::make_pair(pid, fid);

        put(m_dppmb, pid, data_size);
        put(m_rppm, id, ak2k(compute_reference_point(apl0.first, apl1.first, apl2.first)));
        put(m_bpm, id, compute_bbox(apl0.first, apl1.first, apl2.first));

//        std::cout << "Primitive[" << id.first << " " << id.second << "]; "
//                  << "Bbox: [" << get(m_bpm, id) << "] "
//                  << "Point: (" << get(m_rppm, id) << ") "
//                  << "Datum: [" << get(m_dpm, id) << "]" << std::endl;

        Primitive p(id/*, m_dpm, m_rppm*/); // pmaps are external, shared data
        tree.insert(p);
      }
    }

    ++fid;
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_SPLITTING_HELPER_H
