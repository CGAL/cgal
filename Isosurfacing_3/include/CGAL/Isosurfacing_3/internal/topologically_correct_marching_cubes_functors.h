// Copyright (c) 2020 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: ( GPL-3.0-or-later OR LicenseRef-Commercial ) AND MIT
//
// Author(s)     : Julian Stahl
//
// This file incorporates work covered by the following copyright and permission notice:
//
//     MIT License
//
//     Copyright (c) 2020 Roberto Grosso
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy
//     of this software and associated documentation files (the "Software"), to deal
//     in the Software without restriction, including without limitation the rights
//     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//     copies of the Software, and to permit persons to whom the Software is
//     furnished to do so, subject to the following conditions:
//
//     The above copyright notice and this permission notice shall be included in all
//     copies or substantial portions of the Software.
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//     SOFTWARE.
//
//
// The code below uses the version of
// https://github.com/rogrosso/tmc available on 15th of September 2022.
//

#ifndef CGAL_ISOSURFACING_3_INTERNAL_TMC_FUNCTORS_H
#define CGAL_ISOSURFACING_3_INTERNAL_TMC_FUNCTORS_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/marching_cubes_functors.h>
#include <CGAL/Isosurfacing_3/internal/tables.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/concurrent_vector.h>
# include <tbb/concurrent_hash_map.h>
#endif

#include <array>
#include <cmath>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

template <typename Domain_,
          typename PointRange,
          typename PolygonRange>
class TMC_functor
{
private:
  using Domain = Domain_;
  using Point_range = PointRange;
  using Polygon_range = PolygonRange;

  using Geom_traits = typename Domain::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;

  using edge_descriptor = typename Domain::edge_descriptor;
  using cell_descriptor = typename Domain::cell_descriptor;

  using Point_index = std::size_t;
  using Edge_index = std::array<std::size_t, 4>;

#ifdef CGAL_LINKED_WITH_TBB
  struct Hash_compare
  {
    static size_t hash(const Edge_index& key)
    {
      std::size_t res = 17;
      res = res * 31 + std::hash<std::size_t>()(key[0]);
      res = res * 31 + std::hash<std::size_t>()(key[1]);
      res = res * 31 + std::hash<std::size_t>()(key[2]);
      res = res * 31 + std::hash<std::size_t>()(key[3]);
      return res;
    }

    static bool equal(const Edge_index& key1, const Edge_index& key2)
    {
      return key1[0] == key2[0] && key1[1] == key2[1] && key1[2] == key2[2] && key1[3] == key2[3];
    }
  };
#else
  struct Hash
  {
    std::size_t operator()(const Edge_index& key) const
    {
      std::size_t res = 17;
      res = res * 31 + std::hash<std::size_t>()(key[0]);
      res = res * 31 + std::hash<std::size_t>()(key[1]);
      res = res * 31 + std::hash<std::size_t>()(key[2]);
      res = res * 31 + std::hash<std::size_t>()(key[3]);
      return res;
    }
  };
#endif

private:
  const Domain& m_domain;
  FT m_isovalue;

#ifdef CGAL_LINKED_WITH_TBB
  std::atomic<Point_index> m_point_counter;

  tbb::concurrent_vector<Point_3> m_points;
  tbb::concurrent_vector<std::array<Point_index, 3> > m_triangles;

  using Edge_point_map = tbb::concurrent_hash_map<Edge_index, Point_index, Hash_compare>;
  Edge_point_map m_edges;
#else
  Point_index m_point_counter;

  std::vector<Point_3> m_points;
  std::vector<std::array<Point_index, 3> > m_triangles;

  std::unordered_map<Edge_index, Point_index, Hash> m_edges;
#endif

public:
  TMC_functor(const Domain& domain,
              const FT isovalue)
    : m_domain(domain),
      m_isovalue(isovalue),
      m_point_counter(0)
  { }

  void operator()(const cell_descriptor& cell)
  {
    std::array<FT, 8> values;
    std::array<Point_3, 8> corners;
    const std::size_t i_case = get_cell_corners(m_domain, cell, m_isovalue, corners, values);

    // skip empty / full cells
    constexpr std::size_t ones = (1 << 8) - 1;
    if((i_case & ones) == ones || // all bits set
       (i_case & ones) == 0) // no bits set
      return;

    // this is the only difference to the default Marching Cubes
    const int tcm = Cube_table::t_ambig[i_case];
    if(tcm == 105)
    {
      if(p_slice(cell, m_isovalue, corners, values, i_case))
        return;
#ifdef CGAL_ISOSURFACING_3_MC_FUNCTORS_DEBUG
      else
        std::cerr << "WARNING: the result might not be topologically correct" << std::endl;
#endif
    }

    std::array<Point_3, 12> vertices;
    MC_construct_vertices(cell, i_case, corners, values, m_isovalue, m_domain, vertices);

    // @todo improve triangle generation

    // construct triangles
    for(int t=0; t<16; t+=3)
    {
      const int t_index = i_case * 16 + t;

      // if(e_tris_list[t_index] == 0x7f)
      if(Cube_table::triangle_cases[t_index] == -1)
        break;

      // @todo move more of this stuff into the table
      const int eg0 = Cube_table::triangle_cases[t_index + 0];
      const int eg1 = Cube_table::triangle_cases[t_index + 1];
      const int eg2 = Cube_table::triangle_cases[t_index + 2];

      // insert new triangle into list
      const Point_index p0 = add_point(vertices[eg0], compute_edge_index(cell, eg0));
      const Point_index p1 = add_point(vertices[eg1], compute_edge_index(cell, eg1));
      const Point_index p2 = add_point(vertices[eg2], compute_edge_index(cell, eg2));

      add_triangle(p2, p1, p0);
    }
  }

  // returns the created triangle list
  template<typename PR, typename TR>
  void to_triangle_soup(PR& points, TR& triangles) const
  {
    points.insert(points.begin(), m_points.begin(), m_points.end());
    for (const auto& tri : m_triangles) {
      triangles.push_back({ tri[0], tri[1], tri[2] });

      // just a safeguard against arrays of the wrong size
      CGAL_assertion(triangles.back().size() == 3);
    }
  }

private:
  Edge_index compute_edge_index(const cell_descriptor& cell, int edge)
  {
    // edge is in 0 - 11

    // there are 12 edges, assign to each vertex three edges, the global edge numbering
    // consists of 3*global_vertex_id + edge_offset.
    const unsigned long long gei_pattern_ = 670526590282893600ull;

    // the edge global index is given by the vertex global index + the edge offset
    const std::size_t shift = 5 * edge;
    const std::size_t ix = cell[0] + ((gei_pattern_ >> shift) & 1);        // global_edge_id[edge][0];
    const std::size_t iy = cell[1] + ((gei_pattern_ >> (shift + 1)) & 1);  // global_edge_id[edge][1];
    const std::size_t iz = cell[2] + ((gei_pattern_ >> (shift + 2)) & 1);  // global_edge_id[edge][2];
    const std::size_t off_val =      ((gei_pattern_ >> (shift + 3)) & 3);

    return { ix, iy, iz, off_val };
  }

  bool find_point(const Edge_index& e, Point_index& i)
  {
#ifdef CGAL_LINKED_WITH_TBB
    typename Edge_point_map::const_accessor acc;
    if (m_edges.find(acc, e))
    {
      i = acc->second;
      return true;
    }
#else
    auto it = m_edges.find(e);
    if (it != m_edges.end())
    {
      i = it->second;
      return true;
    }
#endif
    return false;
  }

  Point_index add_point(const Point_3& p, const Edge_index& e)
  {

#ifdef CGAL_LINKED_WITH_TBB
    typename Edge_point_map::accessor acc;
    if (!m_edges.insert(acc, e))
      return acc->second;

    const Point_index i = m_point_counter++;
    acc->second = i;
    acc.release();

    m_points.grow_to_at_least(i + 1);
#else
    const Point_index i = m_point_counter;
    auto res = m_edges.insert({e, i});
    if (!res.second)
      return res.first->second;

    ++m_point_counter;
    m_points.resize(i + 1);
#endif

    m_points[i] = p;

    return i;
  }

  Point_index add_point_unchecked(const Point_3& p)
  {
    const Point_index i = m_point_counter++;

#ifdef CGAL_LINKED_WITH_TBB
    m_points.grow_to_at_least(i + 1);
#else
    m_points.resize(i + 1);
#endif
    m_points[i] = p;

    return i;
  }

  void add_triangle(const Point_index p0,
                    const Point_index p1,
                    const Point_index p2)
  {
    m_triangles.push_back({p0, p1, p2});
  }

  bool p_slice(const cell_descriptor& cell,
               const FT i0,
               const std::array<Point_3, 8>& corners,
               const std::array<FT, 8>& values,
               const int i_case)
  {
    typename Geom_traits::Compute_x_3 x_coord = m_domain.geom_traits().compute_x_3_object();
    typename Geom_traits::Compute_y_3 y_coord = m_domain.geom_traits().compute_y_3_object();
    typename Geom_traits::Compute_z_3 z_coord = m_domain.geom_traits().compute_z_3_object();
    typename Geom_traits::Construct_point_3 point = m_domain.geom_traits().construct_point_3_object();

    // code edge end vertices for each of the 12 edges
    const unsigned char l_edges_[12] = {16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98};
    auto get_edge_vertex = [](const int e, unsigned int& v0, unsigned int& v1, const unsigned char l_edges_[12])
    {
      v0 = (unsigned int)(l_edges_[e] & 0xF);
      v1 = (unsigned int)(l_edges_[e] >> 4) & 0xF;
    };

    // A hexahedron has twelve edges, save the intersection of the isosurface with the edge
    // save global edge and global vertex index of isosurface
    std::vector<Point_index> vertices(12);

    // save local coordinate along the edge of intersection point
    std::vector<FT> ecoord(12, FT(0));

    // collect vertices
    unsigned short flag{1};
    for(int eg = 0; eg < 12; ++eg)
    {
      if(flag & Cube_table::intersected_edges[i_case])
      {
        // generate vertex here, do not care at this point if vertex already exists
        unsigned int v0, v1;
        get_edge_vertex(eg, v0, v1, l_edges_);

        // @todo use the domain's interpolation scheme?
        FT l = (i0 - values[v0]) / (values[v1] - values[v0]);
        ecoord[eg] = l;

        // interpolate vertex
        const FT px = (FT(1) - l) * x_coord(corners[v0]) + l * x_coord(corners[v1]);
        const FT py = (FT(1) - l) * y_coord(corners[v0]) + l * y_coord(corners[v1]);
        const FT pz = (FT(1) - l) * z_coord(corners[v0]) + l * z_coord(corners[v1]);

        // add vertex and insert to map
        vertices[eg] = add_point(point(px, py, pz), compute_edge_index(cell, eg));
      }

      // next edge
      flag <<= 1;
    }

    // compute oriented contours
    //
    // A contour consists of segment at the faces connecting the intersection of the
    // isosurface with the edges. For each edge, we store the edge to which the segment
    // is outgoing and the edge from which the segment in coming. Therefore, a contour
    // can be reconstructed by connecting the edges in the direction of the outgoing.
    // The contour is oriented in such a way that the positive vertices are outside.
    // 1. build segments
    // 2. connect segments
    // build up segments
    // set segments map
    unsigned char segm_[12] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
    auto set_segm = [](const int e, const int pos, const int val, unsigned char segm_[12])
    {
      if(pos == 0)
      {
        segm_[e] &= 0xF0;
        segm_[e] |= (unsigned char)val & 0xF;
      }
       else if(pos == 1)
       {
        segm_[e] &= 0xF;
        segm_[e] |= val << 4;
      }
    };

    auto get_segm = [](const int e, const int pos, unsigned char segm_[12]) -> int
    {
      if(pos == 0)
        return int(segm_[e] & 0xF);
      else
        return int((segm_[e] >> 4) & 0xF);
    };

    auto is_segm_set = [](const int e, unsigned char segm_[12]) { return (segm_[e] != 0xFF); };
    auto unset_segm = [](const int e, unsigned char segm_[12]) { segm_[e] = 0xFF; };

    // In order to compute oriented segments, the hexahedron must be flattened.
    // The inside of the faces of the hexahedron must be all on the same
    // side of the flattened hexahedron. This requires changing the order of
    // the edges when reading from the faces
    // code edges at face
    // unsigned short face_e_[6] = { 12816, 30292, 33936, 46754, 34739, 38305 };
    std::array<unsigned short, 6> e_face_{{291, 18277, 18696, 10859, 33719, 38305}};

    // code vertices at face
    // unsigned short face_v_[6] = { 12816, 30292, 21520, 30258, 25632, 30001 };
    std::array<unsigned short, 6> v_face_{{12576, 25717, 5380, 29538, 8292, 30001}};

    // reading edge from face
    auto get_face_e = [e_face_](const int f, const int e) { return ((e_face_[f] >> (4 * e)) & 0xF); };
    auto get_face_v = [v_face_](const int f, const int e) { return ((v_face_[f] >> (4 * e)) & 0xF); };

    // compute oriented segments using the isoline scheme at the faces
    const unsigned int BIT_1 = 1;
    const unsigned int BIT_2 = 2;
    const unsigned int BIT_3 = 4;
    const unsigned int BIT_4 = 8;
    auto asymptotic_decider = [](const FT f0, const FT f1, const FT f2, const FT f3) -> FT
    {
      return (f0 * f3 - f1 * f2) / (f0 + f3 - f1 - f2);
    };

    std::vector<bool> f_flag(6, false);
    for(int f=0; f<6; ++f)
    {
      // classify face
      unsigned int f_case = 0;
      unsigned int v0 = get_face_v(f, 0);
      unsigned int v1 = get_face_v(f, 1);
      unsigned int v2 = get_face_v(f, 2);
      unsigned int v3 = get_face_v(f, 3);
      unsigned int e0 = get_face_e(f, 0);
      unsigned int e1 = get_face_e(f, 1);
      unsigned int e2 = get_face_e(f, 2);
      unsigned int e3 = get_face_e(f, 3);
      FT f0 = values[v0];
      FT f1 = values[v1];
      FT f2 = values[v2];
      FT f3 = values[v3];
      if(f0 >= i0) f_case |= BIT_1;
      if(f1 >= i0) f_case |= BIT_2;
      if(f2 >= i0) f_case |= BIT_3;
      if(f3 >= i0) f_case |= BIT_4;

      switch (f_case)
      {
        case 1:
          set_segm(e0, 0, e3, segm_);
          set_segm(e3, 1, e0, segm_);
        break;
        case 2:
          set_segm(e1, 0, e0, segm_);
          set_segm(e0, 1, e1, segm_);
        break;
        case 3:
          set_segm(e1, 0, e3, segm_);
          set_segm(e3, 1, e1, segm_);
        break;
        case 4:
          set_segm(e3, 0, e2, segm_);
          set_segm(e2, 1, e3, segm_);
        break;
        case 5:
          set_segm(e0, 0, e2, segm_);
          set_segm(e2, 1, e0, segm_);
        break;
        case 6:
        {
          const FT val = asymptotic_decider(f0, f1, f2, f3);
          if(val > i0)
          {
            set_segm(e3, 0, e0, segm_);
            set_segm(e0, 1, e3, segm_);
            set_segm(e1, 0, e2, segm_);
            set_segm(e2, 1, e1, segm_);
          }
          else if(val < i0)
          {
            set_segm(e1, 0, e0, segm_);
            set_segm(e0, 1, e1, segm_);
            set_segm(e3, 0, e2, segm_);
            set_segm(e2, 1, e3, segm_);
          }
          else
          {
            f_flag[f] = true;
            // singular case val == i0, there are no asymptotes
            // check if there is a reasonable triangulation of the face
            const unsigned short e_flag = 0x218;
            const unsigned short bit_1 = 0x1;
            const unsigned short bit_2 = 0x2;
            FT ec0 = ecoord[e0];
            FT ec1 = ecoord[e1];
            FT ec2 = ecoord[e2];
            FT ec3 = ecoord[e3];

            if((e_flag >> (f * 2)) & bit_1)
            {
              ec0 = FT(1) - ec0;
              ec2 = FT(1) - ec2;
            }

            if((e_flag >> (f * 2)) & bit_2)
            {
              ec1 = FT(1) - ec1;
              ec3 = FT(1) - ec3;
            }

            if(ec1 < ec3 && ec0 > ec2)
            {
              set_segm(e1, 0, e0, segm_);
              set_segm(e0, 1, e1, segm_);
              set_segm(e3, 0, e2, segm_);
              set_segm(e2, 1, e3, segm_);
            }
            else if(ec1 > ec3 && ec0 < ec2)
            {
              set_segm(e3, 0, e0, segm_);
              set_segm(e0, 1, e3, segm_);
              set_segm(e1, 0, e2, segm_);
              set_segm(e2, 1, e1, segm_);
            }
            else
            {
              // std::cerr << "ERROR: can't correctly triangulate cell's face\n";
              return false;
            }
          }
        }
        break;
        case 7:
          set_segm(e1, 0, e2, segm_);
          set_segm(e2, 1, e1, segm_);
        break;
        case 8:
          set_segm(e2, 0, e1, segm_);
          set_segm(e1, 1, e2, segm_);
        break;
        case 9:
        {
          const FT val = asymptotic_decider(f0, f1, f2, f3);
          if(val > i0)
          {
            set_segm(e0, 0, e1, segm_);
            set_segm(e1, 1, e0, segm_);
            set_segm(e2, 0, e3, segm_);
            set_segm(e3, 1, e2, segm_);
          }
          else if(val < i0)
          {
            set_segm(e0, 0, e3, segm_);
            set_segm(e3, 1, e0, segm_);
            set_segm(e2, 0, e1, segm_);
            set_segm(e1, 1, e2, segm_);
          }
          else
          {
            f_flag[f] = true;
            // singular case val == i0, there are no asymptotes
            // check if there is a reasonable triangulation of the face
            const unsigned short e_flag = 0x218;
            const unsigned short bit_1 = 0x1;
            const unsigned short bit_2 = 0x2;
            FT ec0 = ecoord[e0];
            FT ec1 = ecoord[e1];
            FT ec2 = ecoord[e2];
            FT ec3 = ecoord[e3];

            if((e_flag >> (f * 2)) & bit_1)
            {
              ec0 = FT(1) - ec0;
              ec2 = FT(1) - ec2;
            }

            if((e_flag >> (f * 2)) & bit_2)
            {
              ec1 = FT(1) - ec1;
              ec3 = FT(1) - ec3;
            }

            if(ec1 < ec3 && ec0 > ec2)
            {
              set_segm(e0, 0, e1, segm_);
              set_segm(e1, 1, e0, segm_);
              set_segm(e2, 0, e3, segm_);
              set_segm(e3, 1, e2, segm_);
            }
            else if(ec1 > ec3 && ec0 < ec2)
            {
              set_segm(e0, 0, e3, segm_);
              set_segm(e3, 1, e0, segm_);
              set_segm(e2, 0, e1, segm_);
              set_segm(e1, 1, e2, segm_);
            }
             else
            {
              // std::cerr << "ERROR: can't correctly triangulate cell's face\n";
              return false;
            }
          }
        }
        break;
        case 10:
          set_segm(e2, 0, e0, segm_);
          set_segm(e0, 1, e2, segm_);
        break;
        case 11:
          set_segm(e2, 0, e3, segm_);
          set_segm(e3, 1, e2, segm_);
        break;
        case 12:
          set_segm(e3, 0, e1, segm_);
          set_segm(e1, 1, e3, segm_);
        break;
        case 13:
          set_segm(e0, 0, e1, segm_);
          set_segm(e1, 1, e0, segm_);
        break;
        case 14:
          set_segm(e3, 0, e0, segm_);
          set_segm(e0, 1, e3, segm_);
        break;
        default:
        break;
      }
    }

    // connect oriented segments into oriented contours
    //
    // closed contours are coded in 64 bit unsigned long long
    // 1) Each entry has 4 bits
    // 2) The first 4 entries are reserved for the size of the contours
    // 3) The next 12 entries are the indices of the edges constituting the contorus
    //    The indices are numbers from 0 to 12
    unsigned long long c_ = 0xFFFFFFFFFFFF0000;

    // in the 4 first bits store size of contours
    auto get_cnt_size = [](const int cnt, unsigned long long& c_) -> size_t
    {
      return size_t((c_ & (0xF << 4 * cnt)) >> 4 * cnt);
    };

    auto set_cnt_size = [](const int cnt, const int size, unsigned long long& c_)
    {
      // unset contour size
      c_ &= ~(0xF << 4 * cnt);
      c_ |= (size << 4 * cnt);
    };

    // set corresponging edge
    auto set_c = [](const int cnt, const int pos, const int val, unsigned long long& c_)
    {
      const unsigned int mask[4] = {0x0, 0xF, 0xFF, 0xFFF};
      const unsigned int c_sz = c_ & mask[cnt];
      const unsigned int e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
      c_ &= ~(((unsigned long long)0xF) << e);
      c_ |= (((unsigned long long)val) << e);
    };

    // read edge from contour
    auto get_c = [](const int cnt, const int pos, unsigned long long c_) -> int
    {
      const unsigned int mask[4] = {0x0, 0xF, 0xFF, 0xFFF};
      const unsigned int c_sz = (unsigned int)(c_ & mask[cnt]);
      const unsigned int e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
      return int((c_ >> e) & 0xF);
    };

    // connect oriented contours
    unsigned int cnt_ = 0;
    for(unsigned int e=0; e<12; ++e)
    {
      if(is_segm_set(e, segm_))
      {
        unsigned int eTo = get_segm(e, 0, segm_);
        unsigned int eIn = get_segm(e, 1, segm_);
        unsigned int eStart = e;
        unsigned int pos = 0;
        set_c(cnt_, pos, eStart, c_);

        while(eTo != eStart)
        {
          pos = pos + 1;
          set_c(cnt_, pos, eTo, c_);
          eIn = eTo;
          eTo = get_segm(eIn, 0, segm_);
          unset_segm(eIn, segm_);
        }

        // set contour length
        set_cnt_size(cnt_, pos + 1, c_);

        // update number of contours
        cnt_ = cnt_ + 1;
      }
    }

    // compute intersection of opposite faces
    //
    // It is sufficient to compute a pair of solutions for one face
    // The other solutions are obtained by evaluating the equations
    // for the common variable
    FT ui[2]{};
    FT vi[2]{};
    FT wi[2]{};
    unsigned char q_sol = 0;
    const FT a = (values[0] - values[1]) * (-values[6] + values[7] + values[4] - values[5]) -
                 (values[4] - values[5]) * (-values[2] + values[3] + values[0] - values[1]);
    const FT b = (i0 - values[0]) * (-values[6] + values[7] + values[4] - values[5]) +
                   (values[0] - values[1]) * (values[6] - values[4]) -
                 (i0 - values[4]) * (-values[2] + values[3] + values[0] - values[1]) -
                   (values[4] - values[5]) * (values[2] - values[0]);
    const FT c = (i0 - values[0]) * (values[6] - values[4]) - (i0 - values[4]) * (values[2] - values[0]);

    FT d = b * b - FT(4) * a * c;
    if(d > 0)
    {
      d = sqrt(d);

      // compute u-coord of solutions
      ui[0] = (-b - d) / (FT(2) * a);
      ui[1] = (-b + d) / (FT(2) * a);

      // compute v-coord of solutions
      FT g1 = values[0] * (FT(1) - ui[0]) + values[1] * ui[0];
      FT g2 = values[2] * (FT(1) - ui[0]) + values[3] * ui[0];
      vi[0] = (i0 - g1) / (g2 - g1);
      if(std::isnan(vi[0]) || std::isinf(vi[0]))
        vi[0] = FT(-1);

      g1 = values[0] * (FT(1) - ui[1]) + values[1] * ui[1];
      g2 = values[2] * (FT(1) - ui[1]) + values[3] * ui[1];
      vi[1] = (i0 - g1) / (g2 - g1);
      if(std::isnan(vi[1]) || std::isinf(vi[1]))
        vi[1] = FT(-1);

      // compute w-coordinates of solutions
      g1 = values[0] * (FT(1) - ui[0]) + values[1] * ui[0];
      g2 = values[4] * (FT(1) - ui[0]) + values[5] * ui[0];
      wi[0] = (i0 - g1) / (g2 - g1);
      if (std::isnan(wi[0]) || std::isinf(wi[0])) wi[0] = FT(-1);
      g1 = values[0] * (FT(1) - ui[1]) + values[1] * ui[1];
      g2 = values[4] * (FT(1) - ui[1]) + values[5] * ui[1];
      wi[1] = (i0 - g1) / (g2 - g1);
      if(std::isnan(wi[1]) || std::isinf(wi[1]))
        wi[1] = FT(-1);

      // correct values for roots of quadratic equations
      // in case the asymptotic decider has failed
      if(f_flag[0]) {  // face 1, w = 0;
        if(wi[0] < wi[1])
          wi[0] = FT(0);
        else
          wi[1] = FT(0);
      }

      if(f_flag[1]) {  // face 2, w = 1
        if(wi[0] > wi[1])
          wi[1] = FT(1);
        else
          wi[1] = FT(1);
      }

      if(f_flag[2]) {  // face 3, v = 0
        if(vi[0] < vi[1])
          vi[0] = FT(0);
        else
          vi[1] = FT(0);
      }

      if(f_flag[3]) {  // face 4, v = 1
        if(vi[0] > vi[1])
          vi[0] = FT(1);
        else
          vi[1] = FT(1);
      }

      if(f_flag[4]) {  // face 5, u = 0
        if(ui[0] < ui[1])
          ui[0] = FT(0);
        else
          ui[1] = FT(0);
      }

      if(f_flag[5]) {  // face 6, u = 1
        if(ui[0] > ui[1])
          ui[0] = FT(1);
        else
          ui[1] = FT(1);
      }

      // check solution intervals
      if(FT(0) < ui[0] && ui[0] < FT(1))
        q_sol |= 1;

      if(0 < ui[1] && ui[1] < 1)
        q_sol |= 2;

      if(0 < vi[0] && vi[0] < 1)
        q_sol |= 4;

      if(0 < vi[1] && vi[1] < 1)
        q_sol |= 8;

      if(0 < wi[0] && wi[0] < 1)
        q_sol |= 16;

      if(0 < wi[1] && wi[1] < 1)
        q_sol |= 32;
    }

    // counts the number of set bits
    auto numberOfSetBits = [](const unsigned char n)
    {
      // C or C++: use uint32_t
      unsigned int b = (unsigned int)(n);
      b = b - ((b >> 1) & 0x55555555);
      b = (b & 0x33333333) + ((b >> 2) & 0x33333333);
      return (((b + (b >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
    };

    // compute the number of solutions to the quadratic equation for a given face
    auto nrQSolFace = [](const unsigned int f, const unsigned char n)
    {
      unsigned int nr = 0;
      switch (f)
      {
        case 0:
          if((n & 0x5) == 0x5) nr = nr + 1;
          if((n & 0xA) == 0xA) nr = nr + 1;
        break;
        case 1:
          if((n & 0x11) == 0x11) nr = nr + 1;
          if((n & 0x22) == 0x22) nr = nr + 1;
        break;
        case 2:
          if((n & 0x18) == 0x18) nr = nr + 1;
          if((n & 0x24) == 0x24) nr = nr + 1;
        break;
      }
      return nr;
    };

    // triangulate contours
    //
    // if all bits are set, then there are three pairs of nontrivial solutions
    // to the quadratic equations. In this case, there is a tunnel or a contour
    // with 12 vertices. If there are three contours, then there is a tunnel and
    // one of the contorus with only three vertices is not part of it.
    if(numberOfSetBits(q_sol) == 6)
    {
      // there are at most three contours
      // Possible cases:
      //  1) a single contour with 12 vertices
      //  2) two contours which build a tunnel
      //  3) three contours, one has only 3 vertices and does not belong to the tunnel

      // construct the six vertices of the inner hexagon
      FT hvt[6][3];
      hvt[0][0] = ui[0];
      hvt[0][1] = vi[0];
      hvt[0][2] = wi[0];
      hvt[1][0] = ui[0];
      hvt[1][1] = vi[0];
      hvt[1][2] = wi[1];
      hvt[2][0] = ui[1];
      hvt[2][1] = vi[0];
      hvt[2][2] = wi[1];
      hvt[3][0] = ui[1];
      hvt[3][1] = vi[1];
      hvt[3][2] = wi[1];
      hvt[4][0] = ui[1];
      hvt[4][1] = vi[1];
      hvt[4][2] = wi[0];
      hvt[5][0] = ui[0];
      hvt[5][1] = vi[1];
      hvt[5][2] = wi[0];

      // construct vertices at intersections with the edges
      auto e_vert = [&ecoord](const int e, const int i) -> FT
      {
        const unsigned int l_coord[3]{1324855, 5299420, 16733440};
        const unsigned char flag = (l_coord[i] >> (2 * e)) & 3;
        if(flag == 3)
          return ecoord[e];
        else
          return FT(flag);
      };

      // if there are three contours, then there is a tunnel and one
      // of the contours is not part of it.
      unsigned char _not_tunnel = 0xF;
      if(cnt_ == 3)
      {
        // loop over the contours
        // triangulate the contour which is not part of
        // the tunnel
        const FT uc_min = (ui[0] < ui[1]) ? ui[0] : ui[1];
        const FT uc_max = (ui[0] < ui[1]) ? ui[1] : ui[0];
        for(int t=0; t < (int)cnt_; ++t)
        {
          if(get_cnt_size(t, c_) == 3)
          {
            FT umin(2);
            FT umax(-2);
            const unsigned int e0 = get_c(t, 0, c_);
            const unsigned int e1 = get_c(t, 1, c_);
            const unsigned int e2 = get_c(t, 2, c_);
            const FT u_e0 = e_vert(e0, 0);
            const FT u_e1 = e_vert(e1, 0);
            const FT u_e2 = e_vert(e2, 0);
            umin = (u_e0 < umin) ? u_e0 : umin;
            umin = (u_e1 < umin) ? u_e1 : umin;
            umin = (u_e2 < umin) ? u_e2 : umin;
            umax = (u_e0 > umax) ? u_e0 : umax;
            umax = (u_e1 > umax) ? u_e1 : umax;
            umax = (u_e2 > umax) ? u_e1 : umax;
            if(uc_min > umax || uc_max < umin)
            {
              // this contour is not part of the tunnel
              _not_tunnel = t;

              add_triangle(vertices[e0], vertices[e1], vertices[e2]);
            }
          }
        }
      }

      // compute vertices of inner hexagon, save new vertices in list and compute and keep
      // global vertices index to build triangle connectivity later on.
      Point_index tg_idx[6];
      for(int i=0; i<6; ++i)
      {
        const FT u = hvt[i][0];
        const FT v = hvt[i][1];
        const FT w = hvt[i][2];
        const FT px = (FT(1) - w) * ((FT(1) - v) * (x_coord(corners[0]) + u * (x_coord(corners[1]) - x_coord(corners[0]))) +
                                               v * (x_coord(corners[2]) + u * (x_coord(corners[3]) - x_coord(corners[2])))) +
                                w * ((FT(1) - v) * (x_coord(corners[4]) + u * (x_coord(corners[5]) - x_coord(corners[4]))) +
                                               v * (x_coord(corners[6]) + u * (x_coord(corners[7]) - x_coord(corners[6]))));
        const FT py = (FT(1) - w) * ((FT(1) - v) * (y_coord(corners[0]) + u * (y_coord(corners[1]) - y_coord(corners[0]))) +
                                               v * (y_coord(corners[2]) + u * (y_coord(corners[3]) - y_coord(corners[2])))) +
                                w * ((FT(1) - v) * (y_coord(corners[4]) + u * (y_coord(corners[5]) - y_coord(corners[4]))) +
                                               v * (y_coord(corners[6]) + u * (y_coord(corners[7]) - y_coord(corners[6]))));
        const FT pz = (FT(1) - w) * ((FT(1) - v) * (z_coord(corners[0]) + u * (z_coord(corners[1]) - z_coord(corners[0]))) +
                                               v * (z_coord(corners[2]) + u * (z_coord(corners[3]) - z_coord(corners[2])))) +
                                w * ((FT(1) - v) * (z_coord(corners[4]) + u * (z_coord(corners[5]) - z_coord(corners[4]))) +
                                               v * (z_coord(corners[6]) + u * (z_coord(corners[7]) - z_coord(corners[6]))));

        tg_idx[i] = add_point_unchecked(point(px, py, pz));
      }

      // triangulate contours with inner hexagon
      unsigned char tcon_[12];
      for(int i=0; i<(int)cnt_; ++i)
      {
        if(_not_tunnel != i)
        {
          // contour belongs to tunnel
          const int cnt_sz = int(get_cnt_size(i, c_));
          for(int r=0; r<cnt_sz; ++r)
          {
            unsigned int index = -1;
            FT dist = std::numeric_limits<FT>::max();
            unsigned int ci = get_c(i, r, c_);
            const FT u_edge = e_vert(ci, 0);
            const FT v_edge = e_vert(ci, 1);
            const FT w_edge = e_vert(ci, 2);
            for(int s=0; s<6; ++s)
            {
              const FT uval = u_edge - hvt[s][0];
              const FT vval = v_edge - hvt[s][1];
              const FT wval = w_edge - hvt[s][2];
              FT val = uval * uval + vval * vval + wval * wval;
              if(dist > val)
              {
                index = s;
                dist = val;
              }
            }

            tcon_[ci] = (unsigned char)(index);
          }

          // correspondence between vertices found
          // create triangles
          // needs some functions
          auto distanceRingIntsModulo = [](const int d1, const int d2)
          {
            const int r = (d1 - d2) < 0 ? d2 - d1 : d1 - d2;
            return (r > 2 ? 6 - r : r);
          };

          auto midpointRingIntModulo = [](const int d1, const int d2)
          {
            const int dmax = (d1 > d2) ? d1 : d2;
            const int dmin = (d1 < d2) ? d1 : d2;
            return ((dmax + 2) % 6 == dmin) ? (dmax + 1) % 6 : (dmax + dmin) / 2;
          };

          for(int r=0; r<cnt_sz; ++r)
          {
            const unsigned int tid1 = get_c(i, r, c_);
            const unsigned int tid2 = get_c(i, ((r + 1) % cnt_sz), c_);
            const unsigned int cid1 = tcon_[tid1];
            const unsigned int cid2 = tcon_[tid2];

            // compute index distance
            const int dst = distanceRingIntsModulo(cid1, cid2);
            switch(dst)
            {
              case 0:
                add_triangle(vertices[tid1], vertices[tid2], tg_idx[cid1]);
              break;
              case 1:
              {
                // measure diagonals
                // triangulate along shortest diagonal
                FT u_edge = e_vert(tid1, 0);
                FT v_edge = e_vert(tid1, 1);
                FT w_edge = e_vert(tid1, 2);
                const FT l1 = (u_edge - hvt[cid2][0]) * (u_edge - hvt[cid2][0]) +
                              (v_edge - hvt[cid2][1]) * (v_edge - hvt[cid2][1]) +
                              (w_edge - hvt[cid2][2]) * (w_edge - hvt[cid2][2]);
                u_edge = e_vert(tid2, 0);
                v_edge = e_vert(tid2, 1);
                w_edge = e_vert(tid2, 2);
                const FT l2 = (u_edge - hvt[cid1][0]) * (u_edge - hvt[cid1][0]) +
                              (v_edge - hvt[cid1][1]) * (v_edge - hvt[cid1][1]) +
                              (w_edge - hvt[cid1][2]) * (w_edge - hvt[cid1][2]);

                if(l1 < l2)
                {
                  add_triangle(vertices[tid1], vertices[tid2], tg_idx[cid2]);
                  add_triangle(vertices[tid1], tg_idx[cid2], tg_idx[cid1]);
                }
                else
                {
                  add_triangle(vertices[tid1], vertices[tid2], tg_idx[cid1]);
                  add_triangle(vertices[tid2], tg_idx[cid2], tg_idx[cid1]);
                }
              }
              break;
              case 2:
              {
                const int cidm = midpointRingIntModulo(cid1, cid2);

                add_triangle(vertices[tid1], vertices[tid2], tg_idx[cidm]);
                add_triangle(vertices[tid1], tg_idx[cidm], tg_idx[cid1]);
                add_triangle(vertices[tid2], tg_idx[cid2], tg_idx[cidm]);
              }
              break;
            } // switch
          } // for loop over the vertices of the contour
        } // if(_not_tunnel)
      } // for loop over contours

      if(cnt_ == 1)
      {
        // there is a single contour
        // triangulate and close inner hexagon
        // triangle must have the correct orientation
        // use asymptotic_decider() to see if positive vertices
        // are separated, in this case orientation must be changed
        const bool s_ = (asymptotic_decider(values[0], values[1], values[2], values[3]) <= i0);
        const bool of_ = (wi[1] < wi[0]) ? s_ : !s_;

        if(!of_)
        {
          add_triangle(tg_idx[0], tg_idx[2], tg_idx[1]);
          add_triangle(tg_idx[2], tg_idx[4], tg_idx[3]);
          add_triangle(tg_idx[0], tg_idx[5], tg_idx[4]);
          add_triangle(tg_idx[0], tg_idx[4], tg_idx[2]);
        }
        else
        {
          add_triangle(tg_idx[0], tg_idx[1], tg_idx[2]);
          add_triangle(tg_idx[2], tg_idx[3], tg_idx[4]);
          add_triangle(tg_idx[0], tg_idx[4], tg_idx[5]);
          add_triangle(tg_idx[0], tg_idx[2], tg_idx[4]);
        }
      }
    }
    else
    {
      // there is no tunnel
      // handle case with no saddle point as simple polygons with 3, 4, 5 or six vertices
      const unsigned char nr_u{(unsigned char)nrQSolFace(0, q_sol)};
      const unsigned char nr_v{(unsigned char)nrQSolFace(1, q_sol)};
      const unsigned char nr_w{(unsigned char)nrQSolFace(2, q_sol)};
      const unsigned char nr_t{(unsigned char)(nr_u + nr_v + nr_w)};
      if(nr_t == nr_u || nr_t == nr_v || nr_t == nr_w)
      {
        // loop over all contours
        for(int i=0; i<(int)cnt_; ++i)
        {
          switch (get_cnt_size(i, c_))
          {
            case 3:
            {
              add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 1, c_)],
                           vertices[get_c(i, 2, c_)]);
            }
            break;
            case 4:
            {
              add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 1, c_)],
                           vertices[get_c(i, 2, c_)]);
              add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 2, c_)],
                           vertices[get_c(i, 3, c_)]);
            }
            break;
            case 5:
            {
              add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 1, c_)],
                           vertices[get_c(i, 2, c_)]);
              add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 2, c_)],
                           vertices[get_c(i, 3, c_)]);
              add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 3, c_)],
                           vertices[get_c(i, 4, c_)]);
            }
            break;
            case 6:
            {
              add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 1, c_)],
                           vertices[get_c(i, 3, c_)]);
              add_triangle(vertices[get_c(i, 1, c_)], vertices[get_c(i, 2, c_)],
                           vertices[get_c(i, 3, c_)]);
              add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 3, c_)],
                           vertices[get_c(i, 4, c_)]);
              add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 4, c_)],
                           vertices[get_c(i, 5, c_)]);
            }
            break;
          } // switch over size of contour
        } // loop over contorus
      }
      else // there are no saddle points
      {
        // there are saddle points
        // fc1 = fs(1, 1)*fs(2, 1) + fs(1, 2)*fs(2, 2);
        // fc2 = fs(1, 1)*fs(3, 1) + fs(1, 2)*fs(3, 2);
        // fc3 = fs(2, 1)*fs(3, 2) + fs(2, 2)*fs(3, 1);
        using uchar = unsigned char;  // @todo

        unsigned char fs[3][2]{{(uchar)(q_sol & 1), (uchar)((q_sol >> 1) & 1)},
                               {(uchar)((q_sol >> 2) & 1), (uchar)((q_sol >> 3) & 1)},
                               {(uchar)((q_sol >> 4) & 1), (uchar)((q_sol >> 5) & 1)}};

        const unsigned char fc1 = fs[0][0] * fs[1][0] + fs[0][1] * fs[1][1];
        const unsigned char fc2 = fs[0][0] * fs[2][0] + fs[0][1] * fs[2][1];
        const unsigned char fc3 = fs[1][0] * fs[2][1] + fs[1][1] * fs[2][0];
        const unsigned char c_faces = fc1 + fc2 + fc3;
        FT ucoord{};
        FT vcoord{};
        FT wcoord{};
        switch(c_faces)
        {
          case 2:
          {
            if(fc1 == 0)
            {
              ucoord = fs[0][0] * ui[0] + fs[0][1] * ui[1];
              vcoord = fs[1][0] * vi[0] + fs[1][1] * vi[1];
              wcoord = fs[1][0] * wi[1] + fs[1][1] * wi[0];
            }
            else if(fc2 == 0)
            {
              ucoord = fs[0][0] * ui[0] + fs[0][1] * ui[1];
              vcoord = fs[0][0] * vi[0] + fs[0][1] * vi[1];
              wcoord = fs[0][0] * wi[1] + fs[0][1] * wi[0];
            }
            else if(fc3 == 0)
            {
              ucoord = fs[1][0] * ui[0] + fs[1][1] * ui[1];
              vcoord = fs[1][0] * vi[0] + fs[1][1] * vi[1];
              wcoord = fs[1][0] * wi[0] + fs[1][1] * wi[1];
            }
          }
          break;
          case 3:
          {
            ucoord = (fs[0][0] * ui[0] + fs[0][1] * ui[1]) / (fs[0][0] + fs[0][1]);
            vcoord = (fs[1][0] * vi[0] + fs[1][1] * vi[1]) / (fs[1][0] + fs[1][1]);
            wcoord = (fs[2][0] * wi[0] + fs[2][1] * wi[1]) / (fs[2][0] + fs[2][1]);
          }
          break;
          case 4:
          {
            const unsigned char nr_u = fs[0][0] + fs[0][1];
            const unsigned char nr_v = fs[1][0] + fs[1][1];
            const unsigned char nr_w = fs[2][0] + fs[2][1];
            if(nr_w == 1)
            {
              ucoord = fs[2][0] * ui[0] + fs[2][1] * ui[1];
              vcoord = fs[2][1] * vi[0] + fs[2][0] * vi[1];
              wcoord = fs[2][0] * wi[0] + fs[2][1] * wi[1];
            }
            else if(nr_v == 1)
            {
              ucoord = fs[1][0] * ui[0] + fs[1][1] * ui[1];
              vcoord = fs[1][0] * vi[0] + fs[1][1] * vi[1];
              wcoord = fs[1][1] * wi[0] + fs[1][0] * wi[1];
            }
            else if(nr_u == 1)
            {
              ucoord = fs[0][0] * ui[0] + fs[0][1] * ui[1];
              vcoord = fs[0][0] * vi[0] + fs[0][1] * vi[1];
              wcoord = fs[0][0] * wi[0] + fs[0][1] * wi[1];
            }
          }
          break;
        } // switch(c_faces)

        // create inner vertex
        const FT px = (FT(1) - wcoord) * ((FT(1) - vcoord) * (x_coord(corners[0]) + ucoord * (x_coord(corners[1]) - x_coord(corners[0]))) +
                                                    vcoord * (x_coord(corners[2]) + ucoord * (x_coord(corners[3]) - x_coord(corners[2])))) +
                                wcoord * ((FT(1) - vcoord) * (x_coord(corners[4]) + ucoord * (x_coord(corners[5]) - x_coord(corners[4]))) +
                                                    vcoord * (x_coord(corners[6]) + ucoord * (x_coord(corners[7]) - x_coord(corners[6]))));
        const FT py = (FT(1) - wcoord) * ((FT(1) - vcoord) * (y_coord(corners[0]) + ucoord * (y_coord(corners[1]) - y_coord(corners[0]))) +
                                                    vcoord * (y_coord(corners[2]) + ucoord * (y_coord(corners[3]) - y_coord(corners[2])))) +
                                wcoord * ((FT(1) - vcoord) * (y_coord(corners[4]) + ucoord * (y_coord(corners[5]) - y_coord(corners[4]))) +
                                                    vcoord * (y_coord(corners[6]) + ucoord * (y_coord(corners[7]) - y_coord(corners[6]))));
        const FT pz = (FT(1) - wcoord) * ((FT(1) - vcoord) * (z_coord(corners[0]) + ucoord * (z_coord(corners[1]) - z_coord(corners[0]))) +
                                                    vcoord * (z_coord(corners[2]) + ucoord * (z_coord(corners[3]) - z_coord(corners[2])))) +
                                wcoord * ((FT(1) - vcoord) * (z_coord(corners[4]) + ucoord * (z_coord(corners[5]) - z_coord(corners[4]))) +
                                                    vcoord * (z_coord(corners[6]) + ucoord * (z_coord(corners[7]) - z_coord(corners[6]))));

        bool pt_used = false;
        Point_index g_index = 0;

        // loop over the contours
        for(int i=0; i<int(cnt_); ++i)
        {
          const unsigned char cnt_sz = (unsigned char)get_cnt_size(i, c_);
          if(cnt_sz == 3)
          {
            add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 1, c_)], vertices[get_c(i, 2, c_)]);
          }
          else
          {
            if (!pt_used)
            {
              pt_used = true;
              g_index = add_point_unchecked(point(px, py, pz));
            }

            for(int t=0; t<cnt_sz; ++t)
            {
              add_triangle(vertices[get_c(i, t, c_)], vertices[get_c(i, (t + 1) % cnt_sz, c_)], g_index);
            }
          }
        }
      }
    }
    return true;
  }
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_TMC_FUNCTORS_H
