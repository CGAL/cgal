// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// Copyright (c) 2010, 2014  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_TRIANGLE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_TRIANGLE_3_DO_INTERSECT_H

// Fast Triangle-Cuboid intersection test, following Tomas Akenine-Moeller description.
// The code looks slightly different from his code because we avoid the translation at
// a minimal cost (and we use C++ ;).

#include <CGAL/Triangle_3.h>
#include <CGAL/Bbox_3.h>

#include <CGAL/Uncertain.h>
#include <CGAL/Intersections_3/internal/Bbox_3_Plane_3_do_intersect.h>

#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K, class Box3>
inline
bool do_bbox_intersect(const typename K::Triangle_3& triangle,
                       const Box3& bbox)
{
  const typename K::Point_3& p = triangle.vertex(0);
  const typename K::Point_3& q = triangle.vertex(1);
  const typename K::Point_3& r = triangle.vertex(2);

  for(int i = 0; i < 3; ++i) {
    if(p[i] <= q[i]) {
      if(q[i] <= r[i]) { // pqr
        if((bbox.max_coord(i) < p[i]) || (bbox.min_coord(i) > r[i]))
          return false;
      }
      else {
        if(p[i] <= r[i]) { // prq
          if(bbox.max_coord(i) < p[i] || bbox.min_coord(i) > q[i])
            return false;
        }
        else { // rpq
          if(bbox.max_coord(i) < r[i] || bbox.min_coord(i) > q[i])
            return false;
        }
      }
    }
    else {
      if(p[i] <= r[i]) { // qpr
        if(bbox.max_coord(i) < q[i] || bbox.min_coord(i) > r[i])
          return false;
      }
      else {
        if(q[i] <= r[i]) { // qrp
          if(bbox.max_coord(i) < q[i] || bbox.min_coord(i) > p[i])
            return false;
        }
        else { // rqp
          if(bbox.max_coord(i) < r[i] || bbox.min_coord(i) > p[i])
            return false;
        }
      }
    }
  }
  return true;
}

// AXE is the axe such that p is orthogonal to it.
// if you do not know it, or if it does not exist,
// use get_min_max without the AXE template parameter
// available in _plane_is_cuboid_do_intersect.h
template <class FT, class Box3, int AXE>
inline
void get_min_max(const FT& px,
                 const FT& py,
                 const FT& pz,
                 const Box3& c,
                 std::array<FT, 3>& p_min,
                 std::array<FT, 3>& p_max)
{
  if(AXE == 0 || px > 0) {
    if(AXE == 1 || py > 0) {
      if(AXE == 2 || pz > 0) {
        p_min = CGAL::make_array<FT>(c.xmin(), c.ymin(),c.zmin());
        p_max = CGAL::make_array<FT>(c.xmax(), c.ymax(),c.zmax());
      }
      else {
        p_min = CGAL::make_array<FT>(c.xmin(), c.ymin(),c.zmax());
        p_max = CGAL::make_array<FT>(c.xmax(), c.ymax(),c.zmin());
      }
    }
    else {
      if(AXE == 2 || pz > 0) {
        p_min = CGAL::make_array<FT>(c.xmin(), c.ymax(),c.zmin());
        p_max = CGAL::make_array<FT>(c.xmax(), c.ymin(),c.zmax());
      }
      else {
        p_min = CGAL::make_array<FT>(c.xmin(), c.ymax(),c.zmax());
        p_max = CGAL::make_array<FT>(c.xmax(), c.ymin(),c.zmin());
      }
    }
  }
  else {
    if(AXE == 1 || py > 0) {
      if(AXE == 2 || pz > 0) {
        p_min = CGAL::make_array<FT>(c.xmax(), c.ymin(),c.zmin());
        p_max = CGAL::make_array<FT>(c.xmin(), c.ymax(),c.zmax());
      }
      else {
        p_min = CGAL::make_array<FT>(c.xmax(), c.ymin(),c.zmax());
        p_max = CGAL::make_array<FT>(c.xmin(), c.ymax(),c.zmin());
      }
    }
    else {
      if(AXE == 2 || pz > 0) {
        p_min = CGAL::make_array<FT>(c.xmax(), c.ymax(),c.zmin());
        p_max = CGAL::make_array<FT>(c.xmin(), c.ymin(),c.zmax());
      }
      else {
        p_min = CGAL::make_array<FT>(c.xmax(), c.ymax(),c.zmax());
        p_max = CGAL::make_array<FT>(c.xmin(), c.ymin(),c.zmin());
      }
    }
  }
}


template <class FT, int SIDE, class Fct>
inline
Uncertain<Sign>
do_axis_intersect_aux_A0(const FT& alpha,
                         const FT& beta,
                         const std::array<std::array<FT,3>, 3>& sides,
                         Fct do_axis_intersect_aux_impl)
{
  return do_axis_intersect_aux_impl(alpha, beta, sides[SIDE][2], sides[SIDE][1]);
}

template <class FT, int SIDE, class Fct>
inline
Uncertain<Sign>
do_axis_intersect_aux_A1(const FT& alpha,
                         const FT& beta,
                         const std::array<std::array<FT,3>, 3>& sides,
                         Fct do_axis_intersect_aux_impl)
{
  return do_axis_intersect_aux_impl(beta, alpha, sides[SIDE][0], sides[SIDE][2]);
}

template <class FT, int SIDE, class Fct>
inline
Uncertain<Sign>
do_axis_intersect_aux_A2(const FT& alpha,
                         const FT& beta,
                         const std::array<std::array<FT,3>, 3>& sides,
                         Fct do_axis_intersect_aux_impl)
{
  return do_axis_intersect_aux_impl(alpha, beta, sides[SIDE][1], sides[SIDE][0]);
}

//given a vector checks whether it is collinear to a base vector
//of the orthonormal frame. return -1 otherwise
template <class FT>
inline
int
collinear_axis(const std::array<FT,3>& side)
{
  if(certainly(side[0] == 0))
  {
    if(certainly(side[1] == 0))
      return 2;
    if(certainly(side[2] == 0))
      return 1;
  }
  else
  {
    if(certainly(side[1] == 0) && certainly(side[2] == 0))
      return 0;
  }
  return -1;
}

template <class FT, class Box3, int AXE, int SIDE, class Fct>
inline
Uncertain<bool> do_axis_intersect(const std::array<std::array<FT, 3>, 3>& triangle,
                                  const std::array<std::array<FT, 3>, 3>& sides,
                                  const Box3& bbox,
                                  Fct do_axis_intersect_aux_impl)
{
  const std::array<FT, 3>& j = triangle[SIDE];
  const std::array<FT, 3>& k = triangle[(SIDE+2)%3];
  const std::array<FT, 3>* j_ptr = &j;
  const std::array<FT, 3>* k_ptr = &k;

  std::array<FT, 3> p_min, p_max;
  get_min_max<FT, Box3, AXE>(AXE==0? 0: AXE==1? sides[SIDE][2]: -sides[SIDE][1],
                             AXE==0? -sides[SIDE][2]: AXE==1? 0: sides[SIDE][0],
                             AXE==0? sides[SIDE][1]: AXE==1? -sides[SIDE][0]:
                             FT(0), bbox, p_min, p_max);

  switch(AXE)
  {
    case 0:
    {
      // t_max >= t_min
      Uncertain<bool> b = do_axis_intersect_aux_A0<FT,SIDE>(k[1]-j[1], k[2]-j[2], sides, do_axis_intersect_aux_impl) != NEGATIVE;
      if(is_indeterminate(b))
        return b;
      if(b)
        std::swap(j_ptr, k_ptr);

      return CGAL_AND((do_axis_intersect_aux_A0<FT,SIDE>(p_min[1]-(*j_ptr)[1], p_min[2]-(*j_ptr)[2],
                                                         sides, do_axis_intersect_aux_impl) != POSITIVE),
                      (do_axis_intersect_aux_A0<FT,SIDE>(p_max[1]-(*k_ptr)[1], p_max[2]-(*k_ptr)[2],
                                                         sides, do_axis_intersect_aux_impl) != NEGATIVE) );
    }
    case 1:
    {
      // t_max >= t_min
      Uncertain<bool> b = do_axis_intersect_aux_A1<FT,SIDE>(k[0]-j[0], k[2]-j[2], sides, do_axis_intersect_aux_impl) != NEGATIVE;
      if(is_indeterminate(b))
        return b;
      if(b)
        std::swap(j_ptr, k_ptr);

      return CGAL_AND((do_axis_intersect_aux_A1<FT,SIDE>(p_min[0]-(*j_ptr)[0], p_min[2]-(*j_ptr)[2],
                                                         sides, do_axis_intersect_aux_impl) != POSITIVE),
                      (do_axis_intersect_aux_A1<FT,SIDE>(p_max[0]-(*k_ptr)[0], p_max[2]-(*k_ptr)[2],
                                                         sides, do_axis_intersect_aux_impl) != NEGATIVE) );
    }
    case 2:
    {
      // t_max >= t_min
      Uncertain<bool> b = do_axis_intersect_aux_A2<FT,SIDE>(k[0]-j[0], k[1]-j[1], sides, do_axis_intersect_aux_impl) != NEGATIVE;
      if(is_indeterminate(b))
        return b;
      if(b)
        std::swap(j_ptr, k_ptr);

      return CGAL_AND((do_axis_intersect_aux_A2<FT,SIDE>(p_min[0]-(*j_ptr)[0], p_min[1]-(*j_ptr)[1],
                                                         sides, do_axis_intersect_aux_impl) != POSITIVE),
                      (do_axis_intersect_aux_A2<FT,SIDE>(p_max[0]-(*k_ptr)[0], p_max[1]-(*k_ptr)[1],
                                                         sides, do_axis_intersect_aux_impl) != NEGATIVE) );
    }
    default: // Should not happen
      CGAL_error();
      return make_uncertain(false);
  }
}

template <class FT, class Box3, class Fct>
Uncertain<bool>
do_intersect_bbox_or_iso_cuboid_impl(const std::array< std::array<FT, 3>, 3>& triangle,
                                     const Box3& bbox,
                                     Fct do_axis_intersect_aux_impl)
{
  std::array< std::array<FT, 3>, 3> sides =
  {{
    {triangle[1][0] - triangle[0][0], triangle[1][1] - triangle[0][1], triangle[1][2] - triangle[0][2]},
    {triangle[2][0] - triangle[1][0], triangle[2][1] - triangle[1][1], triangle[2][2] - triangle[1][2]},
    {triangle[0][0] - triangle[2][0], triangle[0][1] - triangle[2][1], triangle[0][2] - triangle[2][2]}
  }};

  int forbidden_axis = -1;
  int forbidden_size = -1;
  //determine whether one vector is collinear with an axis
  int tmp = collinear_axis<FT>(sides[0]);
  if(tmp != -1)
  {
    forbidden_axis = tmp;
    forbidden_size = 0;
  }
  else
  {
    tmp = collinear_axis<FT>(sides[1]);
    if(tmp != -1)
    {
      forbidden_axis = tmp;
      forbidden_size = 1;
    }
    else
    {
      tmp = collinear_axis<FT>(sides[2]);
      if(tmp != -1)
      {
        forbidden_axis = tmp;
        forbidden_size = 2;
      }
    }
  }

  // Create a "certainly true"
  Uncertain<bool> ind_or_true = make_uncertain(true);

  if(forbidden_axis != 0)
  {
    if(forbidden_size != 0)
    {
      Uncertain<bool> b = do_axis_intersect<FT,Box3,0,0>(triangle, sides, bbox, do_axis_intersect_aux_impl);
      if(is_indeterminate(b))
        ind_or_true = b;
       else if(!b)
        return false;
    }

    if(forbidden_size != 1)
    {
      Uncertain<bool> b = do_axis_intersect<FT,Box3,0,1>(triangle, sides, bbox, do_axis_intersect_aux_impl);
      if(is_indeterminate(b))
        ind_or_true = b;
       else if(!b)
        return false;
    }

    if(forbidden_size != 2)
    {
      Uncertain<bool> b = do_axis_intersect<FT,Box3,0,2>(triangle, sides, bbox, do_axis_intersect_aux_impl);
      if(is_indeterminate(b))
        ind_or_true = b;
       else if(!b)
        return false;
    }
  }

  if(forbidden_axis != 1)
  {
    if(forbidden_size != 0)
    {
      Uncertain<bool> b = do_axis_intersect<FT,Box3,1,0>(triangle, sides, bbox, do_axis_intersect_aux_impl);
      if(is_indeterminate(b))
        ind_or_true = b;
       else if(!b)
        return false;
    }

    if(forbidden_size != 1)
    {
      Uncertain<bool> b = do_axis_intersect<FT,Box3,1,1>(triangle, sides, bbox, do_axis_intersect_aux_impl);
      if(is_indeterminate(b))
        ind_or_true = b;
       else if(!b)
        return false;
    }

    if(forbidden_size != 2)
    {
      Uncertain<bool> b = do_axis_intersect<FT,Box3,1,2>(triangle, sides, bbox, do_axis_intersect_aux_impl);
      if(is_indeterminate(b))
        ind_or_true = b;
       else if(!b)
        return false;
    }
  }

  if(forbidden_axis != 2)
  {
    if(forbidden_size != 0)
    {
      Uncertain<bool> b = do_axis_intersect<FT,Box3,2,0>(triangle, sides, bbox, do_axis_intersect_aux_impl);
      if(is_indeterminate(b))
        ind_or_true = b;
       else if(!b)
        return false;
    }

    if(forbidden_size != 1)
    {
      Uncertain<bool> b = do_axis_intersect<FT,Box3,2,1>(triangle, sides, bbox, do_axis_intersect_aux_impl);
      if(is_indeterminate(b))
        ind_or_true = b;
      else if(!b)
        return false;
    }

    if(forbidden_size != 2)
    {
      Uncertain<bool> b = do_axis_intersect<FT,Box3,2,2>(triangle, sides, bbox, do_axis_intersect_aux_impl);
      if(is_indeterminate(b))
        ind_or_true = b;
      else if(!b)
        return false;
    }
  }

  return ind_or_true;
}

template <class K, class Box3>
bool do_intersect_bbox_or_iso_cuboid(const typename K::Triangle_3& a_triangle,
                                     const Box3& a_bbox,
                                     const K& k)
{
  if(certainly_not(do_bbox_intersect<K>(a_triangle, a_bbox)))
    return false;

  if(certainly_not(do_intersect(a_triangle.supporting_plane(), a_bbox, k)))
    return false;

  typedef typename K::FT FT;
  auto do_axis_intersect_aux_impl = [](const FT& alpha,
                                       const FT& beta,
                                       const FT& c_alpha,
                                       const FT& c_beta) -> Uncertain<Sign>
  {
    return CGAL::sign(- c_alpha * alpha + c_beta * beta);
  };

  std::array< std::array<FT, 3>, 3> triangle =
  {{
    { a_triangle[0][0], a_triangle[0][1], a_triangle[0][2] },
    { a_triangle[1][0], a_triangle[1][1], a_triangle[1][2] },
    { a_triangle[2][0], a_triangle[2][1], a_triangle[2][2] }
  }};

  // exception will be thrown in case the output is indeterminate
  return do_intersect_bbox_or_iso_cuboid_impl<FT>(triangle, a_bbox, do_axis_intersect_aux_impl);
}

template <class K>
bool do_intersect(const typename K::Triangle_3& triangle,
                  const CGAL::Bbox_3& bbox,
                  const K& k)
{
  return do_intersect_bbox_or_iso_cuboid(triangle, bbox, k);
}

template <class K>
bool do_intersect(const CGAL::Bbox_3& bbox,
                  const typename K::Triangle_3& triangle,
                  const K& k)
{
  return do_intersect_bbox_or_iso_cuboid(triangle, bbox, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_TRIANGLE_3_DO_INTERSECT_H
