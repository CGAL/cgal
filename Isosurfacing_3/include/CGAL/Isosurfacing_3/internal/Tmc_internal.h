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

#ifndef CGAL_TMC_INTERNAL_TMC_H
#define CGAL_TMC_INTERNAL_TMC_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/Marching_cubes_3_internal.h>
#include <CGAL/Isosurfacing_3/internal/Tables.h>

#include <cmath>
#include <array>
#include <atomic>
#include <map>
#include <mutex>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

template <class Domain_, class PointRange, class PolygonRange>
class TMC_functor {
private:
    typedef Domain_ Domain;
    typedef PointRange Point_range;
    typedef PolygonRange Polygon_range;

    typedef typename Domain::FT FT;
    typedef typename Domain::Point Point;
    typedef typename Domain::Vector Vector;
    typedef typename Domain::Edge_descriptor Edge_descriptor;
    typedef typename Domain::Cell_descriptor Cell_descriptor;

    typedef unsigned int uint;

public:
    TMC_functor(const Domain& domain, const FT isovalue, Point_range& points, Polygon_range& polygons)
        : domain(domain), isovalue(isovalue), points(points), polygons(polygons) {}

    void operator()(const Cell_descriptor& cell) {

        FT values[8];
        Point corners[8];
        const int i_case = get_cell_corners(domain, cell, isovalue, corners, values);

        const int all_bits_set = (1 << (8 + 1)) - 1;  // last 8 bits are 1
        if (Cube_table::intersected_edges[i_case] == 0 || Cube_table::intersected_edges[i_case] == all_bits_set) {
            return;
        }

        // this is the only difference to mc
        int tcm = (int)Cube_table::t_ambig[i_case];
        if (tcm == 105) {
            p_slice(cell, isovalue, values, corners, i_case);
            return;
        }

        std::array<Point, 12> vertices;
        mc_construct_vertices(domain.cell_edges(cell), isovalue, i_case, corners, values, vertices);

        // TODO: improve triangle generation
        // construct triangles
        std::lock_guard<std::mutex> lock(mutex);
        for (int t = 0; t < 16; t += 3) {

            const int t_index = i_case * 16 + t;
            // if (e_tris_list[t_index] == 0x7f)
            if (Cube_table::triangle_cases[t_index] == -1) break;

            const int eg0 = Cube_table::triangle_cases[t_index + 0];  // TODO: move more of this stuff into the table
            const int eg1 = Cube_table::triangle_cases[t_index + 1];
            const int eg2 = Cube_table::triangle_cases[t_index + 2];

            const std::size_t p0_idx = points.size();

            points.push_back(vertices[eg0]);
            points.push_back(vertices[eg1]);
            points.push_back(vertices[eg2]);

            // insert new triangle in list
            polygons.push_back({});
            auto& triangle = polygons.back();

            triangle.push_back(p0_idx + 2);
            triangle.push_back(p0_idx + 1);
            triangle.push_back(p0_idx + 0);
        }
    }

    void add_triangle(const std::size_t p0, const std::size_t p1, const std::size_t p2) {
        std::lock_guard<std::mutex> lock(mutex);

        polygons.push_back({});
        auto& triangle = polygons.back();

        triangle.push_back(p0);
        triangle.push_back(p1);
        triangle.push_back(p2);
    }

    void p_slice(const Cell_descriptor& cell, const double i0, FT* values, Point* corners, const int i_case) {
        // there are 12 edges, assign to each vertex three edges, the global edge numbering
        // consist of 3*global_vertex_id + edge_offset.
        const unsigned long long gei_pattern_ = 670526590282893600ull;

        // code edge end vertices for each of the 12 edges
        const unsigned char l_edges_[12] = {16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98};
        auto get_edge_vertex = [](const int e, unsigned int& v0, unsigned int& v1, const unsigned char l_edges_[12]) {
            v0 = (unsigned int)(l_edges_[e] & 0xF);
            v1 = (unsigned int)(l_edges_[e] >> 4) & 0xF;
        };

        // A hexahedron has twelve edges, save the intersection of the isosurface with the edge
        // save global edge and global vertex index of isosurface
        std::vector<std::size_t> vertices(12);
        // save loca coordinate along the edge of intersection point
        std::vector<FT> ecoord{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

        // collect vertices
        unsigned short flag{1};
        for (int eg = 0; eg < 12; eg++) {
            if (flag & Cube_table::intersected_edges[i_case]) {
                // the edge global index is given by the vertex global index + the edge offset
                // uint shift = 5 * eg;
                // const int ix = i_index + (int)((gei_pattern_ >> shift) & 1);        // global_edge_id[eg][0];
                // const int iy = j_index + (int)((gei_pattern_ >> (shift + 1)) & 1);  // global_edge_id[eg][1];
                // const int iz = k_index + (int)((gei_pattern_ >> (shift + 2)) & 1);  // global_edge_id[eg][2];
                // const int off_val = (int)((gei_pattern_ >> (shift + 3)) & 3);

                // int g_edg = int(m_cell_shift_factor * m_ugrid.global_index(ix, iy, iz) + off_val);

                // generate vertex here, do not care at this point if vertex already exist
                uint v0, v1;
                get_edge_vertex(eg, v0, v1, l_edges_);

                double l = (i0 - values[v0]) / (values[v1] - values[v0]);
                ecoord[eg] = l;
                // interpolate vertex
                const FT px = (1 - l) * corners[v0][0] + l * corners[v1][0];
                const FT py = (1 - l) * corners[v0][1] + l * corners[v1][1];
                const FT pz = (1 - l) * corners[v0][2] + l * corners[v1][2];

                // set vertex in map
                // set vertex index
                // auto s_index = m_vertices.find(vertices[eg].g_edg);
                // if (s_index == m_vertices.end()) {
                const int g_idx = (int)points.size();
                vertices[eg] = g_idx;
                //    m_vertices[vertices[eg].g_edg] = g_idx;
                points.push_back(Point(px, py, pz));
                //} else {
                //    vertices[eg] = s_index->second;
                //}
            }
            /*else {
            e_set[eg] = false;
            }*/
            // next edge
            flag <<= 1;
        }

        // compute oriented contours
        // A countour consists of segment at the faces connecting the intersection of the
        // isosurface with the edges. For each edge we store the edge to which the segment
        // is outgoing and the edge from which the segment in comming. Therefore a contour
        // cab be reconstructed by connecting the edges in the direccion of the outgoing.
        // The contour is oriented in such a way, that the positive vertices are outside.
        // 1. build segments
        // 2. connect segments
        // build up segments
        // set segments map
        unsigned char segm_[12] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
        auto set_segm = [](const int e, const int pos, const int val, unsigned char segm_[12]) {
            if (pos == 0) {
                segm_[e] &= 0xF0;
                segm_[e] |= (unsigned char)val & 0xF;
            } else if (pos == 1) {
                segm_[e] &= 0xF;
                segm_[e] |= val << 4;
            }
        };
        auto get_segm = [](const int e, const int pos, unsigned char segm_[12]) {
            if (pos == 0)
                return (int)(segm_[e] & 0xF);
            else
                return (int)((segm_[e] >> 4) & 0xF);
        };
        auto is_segm_set = [](const int e, unsigned char segm_[12]) { return (segm_[e] != 0xFF); };
        auto unset_segm = [](const int e, unsigned char segm_[12]) { segm_[e] = 0xFF; };
        // In order to compute oriented segments, the hexahedron has to be flatten.
        // The insides of the faces of the hexahedron have to be all at the same
        // side of the flattend hexa. This requires changing the order of the
        // edges when reading from the faces
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
        auto asymptotic_decider = [](const double f0, const double f1, const double f2, const double f3) {
            return (f0 * f3 - f1 * f2) / (f0 + f3 - f1 - f2);
        };
        std::vector<bool> f_flag(6, false);
        for (int f = 0; f < 6; f++) {
            // classify face
            unsigned int f_case{0};
            uint v0 = get_face_v(f, 0);
            uint v1 = get_face_v(f, 1);
            uint v2 = get_face_v(f, 2);
            uint v3 = get_face_v(f, 3);
            uint e0 = get_face_e(f, 0);
            uint e1 = get_face_e(f, 1);
            uint e2 = get_face_e(f, 2);
            uint e3 = get_face_e(f, 3);
            double f0 = values[v0];
            double f1 = values[v1];
            double f2 = values[v2];
            double f3 = values[v3];
            if (f0 >= i0) f_case |= BIT_1;
            if (f1 >= i0) f_case |= BIT_2;
            if (f2 >= i0) f_case |= BIT_3;
            if (f3 >= i0) f_case |= BIT_4;
            switch (f_case) {
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
                case 6: {
                    const double val = asymptotic_decider(f0, f1, f2, f3);
                    if (val > i0) {
                        set_segm(e3, 0, e0, segm_);
                        set_segm(e0, 1, e3, segm_);
                        set_segm(e1, 0, e2, segm_);
                        set_segm(e2, 1, e1, segm_);
                    } else if (val < i0) {
                        set_segm(e1, 0, e0, segm_);
                        set_segm(e0, 1, e1, segm_);
                        set_segm(e3, 0, e2, segm_);
                        set_segm(e2, 1, e3, segm_);
                    } else {
                        f_flag[f] = true;
                        // singular case val == i0, there are no asymptotes
                        // check if there is a reasonable triangulation of the face
                        unsigned short e_flag = 0x218;
                        unsigned short bit_1 = 0x1;
                        unsigned short bit_2 = 0x2;
                        double ec0 = ecoord[e0];
                        double ec1 = ecoord[e1];
                        double ec2 = ecoord[e2];
                        double ec3 = ecoord[e3];
                        if ((e_flag >> (f * 2)) & bit_1) {
                            ec0 = 1 - ec0;
                            ec2 = 1 - ec2;
                        }
                        if ((e_flag >> (f * 2)) & bit_2) {
                            ec1 = 1 - ec1;
                            ec3 = 1 - ec3;
                        }
                        if (ec1 < ec3 && ec0 > ec2) {
                            set_segm(e1, 0, e0, segm_);
                            set_segm(e0, 1, e1, segm_);
                            set_segm(e3, 0, e2, segm_);
                            set_segm(e2, 1, e3, segm_);
                        } else if (ec1 > ec3 && ec0 < ec2) {
                            set_segm(e3, 0, e0, segm_);
                            set_segm(e0, 1, e3, segm_);
                            set_segm(e1, 0, e2, segm_);
                            set_segm(e2, 1, e1, segm_);
                        } else {
                            std::cerr << "ERROR: can't correctly triangulate cell's face\n";
                            return;
                        }
                    }
                } break;
                case 7:
                    set_segm(e1, 0, e2, segm_);
                    set_segm(e2, 1, e1, segm_);
                    break;
                case 8:
                    set_segm(e2, 0, e1, segm_);
                    set_segm(e1, 1, e2, segm_);
                    break;
                case 9: {
                    const double val = asymptotic_decider(f0, f1, f2, f3);
                    if (val > i0) {
                        set_segm(e0, 0, e1, segm_);
                        set_segm(e1, 1, e0, segm_);
                        set_segm(e2, 0, e3, segm_);
                        set_segm(e3, 1, e2, segm_);
                    } else if (val < i0) {
                        set_segm(e0, 0, e3, segm_);
                        set_segm(e3, 1, e0, segm_);
                        set_segm(e2, 0, e1, segm_);
                        set_segm(e1, 1, e2, segm_);
                    } else {
                        f_flag[f] = true;
                        // singular case val == i0, there are no asymptotes
                        // check if there is a reasonable triangulation of the face
                        unsigned short e_flag = 0x218;
                        unsigned short bit_1 = 0x1;
                        unsigned short bit_2 = 0x2;
                        double ec0 = ecoord[e0];
                        double ec1 = ecoord[e1];
                        double ec2 = ecoord[e2];
                        double ec3 = ecoord[e3];
                        if ((e_flag >> (f * 2)) & bit_1) {
                            ec0 = 1 - ec0;
                            ec2 = 1 - ec2;
                        }
                        if ((e_flag >> (f * 2)) & bit_2) {
                            ec1 = 1 - ec1;
                            ec3 = 1 - ec3;
                        }
                        if (ec1 < ec3 && ec0 > ec2) {
                            set_segm(e0, 0, e1, segm_);
                            set_segm(e1, 1, e0, segm_);
                            set_segm(e2, 0, e3, segm_);
                            set_segm(e3, 1, e2, segm_);
                        } else if (ec1 > ec3 && ec0 < ec2) {
                            set_segm(e0, 0, e3, segm_);
                            set_segm(e3, 1, e0, segm_);
                            set_segm(e2, 0, e1, segm_);
                            set_segm(e1, 1, e2, segm_);
                        } else {
                            std::cerr << "ERROR: can't correctly triangulate cell's face\n";
                            return;
                        }
                    }
                } break;
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
        // closed contours are coded in 64 bit unsigned long long
        // 1) Each entry has 4 bits
        // 2) The first 4 entries are reserved for the size of the contours
        // 3) The next 12 entries are the indices of the edges constituting the contorus
        //    The indices are numbers from 0 to 12
        unsigned long long c_ = 0xFFFFFFFFFFFF0000;
        // in the 4 first bits store size of contours
        auto get_cnt_size = [](const int cnt, unsigned long long& c_) {
            return (size_t)((c_ & (0xF << 4 * cnt)) >> 4 * cnt);
        };
        auto set_cnt_size = [](const int cnt, const int size, unsigned long long& c_) {
            // unset contour size
            c_ &= ~(0xF << 4 * cnt);
            c_ |= (size << 4 * cnt);
        };
        // set corresponging edge
        auto set_c = [](const int cnt, const int pos, const int val, unsigned long long& c_) {
            const uint mask[4] = {0x0, 0xF, 0xFF, 0xFFF};
            const uint c_sz = c_ & mask[cnt];
            const uint e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
            c_ &= ~(((unsigned long long)0xF) << e);
            c_ |= (((unsigned long long)val) << e);
        };
        // read edge from contour
        auto get_c = [](const int cnt, const int pos, unsigned long long c_) {
            const uint mask[4] = {0x0, 0xF, 0xFF, 0xFFF};
            const uint c_sz = (uint)(c_ & mask[cnt]);
            const uint e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
            return (int)((c_ >> e) & 0xF);
        };


        // connect oriented contours
        uint cnt_{0};
        for (uint e = 0; e < 12; e++) {
            if (is_segm_set(e, segm_)) {
                uint eTo = get_segm(e, 0, segm_);
                uint eIn = get_segm(e, 1, segm_);
                uint eStart = e;
                uint pos = 0;
                set_c(cnt_, pos, eStart, c_);
                while (eTo != eStart) {
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
        // It is enough to compute a pair of solutions for one face
        // The other solutions are obtained by evaluating the equations
        // for the common variable
        double ui[2]{};
        double vi[2]{};
        double wi[2]{};
        unsigned char q_sol{0};
        const double a = (values[0] - values[1]) * (-values[6] + values[7] + values[4] - values[5]) -
                         (values[4] - values[5]) * (-values[2] + values[3] + values[0] - values[1]);
        const double b = (i0 - values[0]) * (-values[6] + values[7] + values[4] - values[5]) +
                         (values[0] - values[1]) * (values[6] - values[4]) -
                         (i0 - values[4]) * (-values[2] + values[3] + values[0] - values[1]) -
                         (values[4] - values[5]) * (values[2] - values[0]);
        const double c = (i0 - values[0]) * (values[6] - values[4]) - (i0 - values[4]) * (values[2] - values[0]);
        ;
        double d = b * b - 4 * a * c;
        if (d > 0) {
            d = std::sqrt(d);
            // compute u-coord of solutions
            ui[0] = (-b - d) / (2 * a);
            ui[1] = (-b + d) / (2 * a);
            // compute v-coord of solutions
            double g1 = values[0] * (1 - ui[0]) + values[1] * ui[0];
            double g2 = values[2] * (1 - ui[0]) + values[3] * ui[0];
            vi[0] = (i0 - g1) / (g2 - g1);
            if (std::isnan(vi[0]) || std::isinf(vi[0])) vi[0] = -1.f;
            g1 = values[0] * (1 - ui[1]) + values[1] * ui[1];
            g2 = values[2] * (1 - ui[1]) + values[3] * ui[1];
            vi[1] = (i0 - g1) / (g2 - g1);
            if (std::isnan(vi[1]) || std::isinf(vi[1])) vi[1] = -1.f;
            // compute w-coordinates of solutions
            g1 = values[0] * (1 - ui[0]) + values[1] * ui[0];
            g2 = values[4] * (1 - ui[0]) + values[5] * ui[0];
            wi[0] = (i0 - g1) / (g2 - g1);
            if (std::isnan(wi[0]) || std::isinf(wi[0])) wi[0] = -1.f;
            g1 = values[0] * (1 - ui[1]) + values[1] * ui[1];
            g2 = values[4] * (1 - ui[1]) + values[5] * ui[1];
            wi[1] = (i0 - g1) / (g2 - g1);
            if (std::isnan(wi[1]) || std::isinf(wi[1])) wi[1] = -1.f;
            // correct values for roots of quadratic equations
            // in case the asymptotic decider has failed
            if (f_flag[0] == true) {  // face 1, w = 0;
                if (wi[0] < wi[1])
                    wi[0] = 0;
                else
                    wi[1] = 0;
            }
            if (f_flag[1] == true) {  // face 2, w = 1
                if (wi[0] > wi[1])
                    wi[1] = 1;
                else
                    wi[1] = 1;
            }
            if (f_flag[2] == true) {  // face 3, v = 0
                if (vi[0] < vi[1])
                    vi[0] = 0;
                else
                    vi[1] = 0;
            }
            if (f_flag[3] == true) {  // face 4, v = 1
                if (vi[0] > vi[1])
                    vi[0] = 1;
                else
                    vi[1] = 1;
            }
            if (f_flag[4] == true) {  // face 5, u = 0
                if (ui[0] < ui[1])
                    ui[0] = 0;
                else
                    ui[1] = 0;
            }
            if (f_flag[5] == true) {  // face 6, u = 1
                if (ui[0] > ui[1])
                    ui[0] = 1;
                else
                    ui[1] = 1;
            }

            // check solution intervals
            if (0 < ui[0] && ui[0] < 1) {
                q_sol |= 1;
            }
            if (0 < ui[1] && ui[1] < 1) {
                q_sol |= 2;
            }
            if (0 < vi[0] && vi[0] < 1) {
                q_sol |= 4;
            }
            if (0 < vi[1] && vi[1] < 1) {
                q_sol |= 8;
            }
            if (0 < wi[0] && wi[0] < 1) {
                q_sol |= 16;
            }
            if (0 < wi[1] && wi[1] < 1) {
                q_sol |= 32;
            }
        }

        //
        // count the number of set bits
        auto numberOfSetBits = [](const unsigned char n) {
            // C or C++: use uint32_t
            uint b = (uint)n;
            b = b - ((b >> 1) & 0x55555555);
            b = (b & 0x33333333) + ((b >> 2) & 0x33333333);
            return (((b + (b >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
        };
        // compute the number of solutions to the quadratic equation for a given face
        auto nrQSolFace = [](const uint f, const unsigned char n) {
            uint nr{0};
            switch (f) {
                case 0:
                    if ((n & 0x5) == 0x5) nr = nr + 1;
                    if ((n & 0xA) == 0xA) nr = nr + 1;
                    break;
                case 1:
                    if ((n & 0x11) == 0x11) nr = nr + 1;
                    if ((n & 0x22) == 0x22) nr = nr + 1;
                    break;
                case 2:
                    if ((n & 0x18) == 0x18) nr = nr + 1;
                    if ((n & 0x24) == 0x24) nr = nr + 1;
                    break;
            }
            return nr;
        };


        // triangulate contours
        // if all bits are set, then there are three pairs of nontrivial solutions
        // to the quadratic equations. In this case, there is a tunnel or a contour
        // with 12 vertices. If there are three contours, then there is a tunnel and
        // one of the contorus with only three vertices is not part of it.
        if (numberOfSetBits(q_sol) == 6) {
            // there are at most three contours
            // Possible cases:
            //  1) a single contour with 12 vertices
            //  2) two contours which build a tunnel
            //  3) three contours, one has only 3 vertices and does not belong to the tunnel

            // construct the six vertices of the inner hexagon
            double hvt[6][3];
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
            auto e_vert = [&ecoord](const int e, const int i) {
                const unsigned int l_coord[3]{1324855, 5299420, 16733440};
                unsigned char flag = (l_coord[i] >> (2 * e)) & 3;
                if (flag == 3)
                    return ecoord[e];
                else
                    return (FT)(flag);
            };

            // if there are three contours, then there is a tunnel and one
            // of the contours is not part of it.
            unsigned char _not_tunnel = 0xF;
            if (cnt_ == 3) {
                // loop over the contorus
                // triangulate the contour which is not part of
                // the tunnel
                const double uc_min = (ui[0] < ui[1]) ? ui[0] : ui[1];
                const double uc_max = (ui[0] < ui[1]) ? ui[1] : ui[0];
                for (int t = 0; t < (int)cnt_; t++) {
                    if (get_cnt_size(t, c_) == 3) {
                        double umin = 2;
                        double umax = -2;
                        uint e0 = get_c(t, 0, c_);
                        uint e1 = get_c(t, 1, c_);
                        uint e2 = get_c(t, 2, c_);
                        const double u_e0 = e_vert(e0, 0);
                        const double u_e1 = e_vert(e1, 0);
                        const double u_e2 = e_vert(e2, 0);
                        umin = (u_e0 < umin) ? u_e0 : umin;
                        umin = (u_e1 < umin) ? u_e1 : umin;
                        umin = (u_e2 < umin) ? u_e2 : umin;
                        umax = (u_e0 > umax) ? u_e0 : umax;
                        umax = (u_e1 > umax) ? u_e1 : umax;
                        umax = (u_e2 > umax) ? u_e1 : umax;
                        if (uc_min > umax || uc_max < umin) {
                            // this contour is not part of the tunnel
                            _not_tunnel = t;

                            add_triangle(vertices[e0], vertices[e1], vertices[e2]);
                        }
                    }
                }
            }

            // compute vertices of inner hexagon, save new vertices in list and compute and keep
            // global vertice index to build triangle connectivity later on.
            uint tg_idx[6];
            for (int i = 0; i < 6; i++) {

                const double u = hvt[i][0];
                const double v = hvt[i][1];
                const double w = hvt[i][2];
                const FT px = (1 - w) * ((1 - v) * (corners[0][0] + u * (corners[1][0] - corners[0][0])) +
                                         v * (corners[2][0] + u * (corners[3][0] - corners[2][0]))) +
                              w * ((1 - v) * (corners[4][0] + u * (corners[5][0] - corners[4][0])) +
                                   v * (corners[6][0] + u * (corners[7][0] - corners[6][0])));
                const FT py = (1 - w) * ((1 - v) * (corners[0][1] + u * (corners[1][1] - corners[0][1])) +
                                         v * (corners[2][1] + u * (corners[3][1] - corners[2][1]))) +
                              w * ((1 - v) * (corners[4][1] + u * (corners[5][1] - corners[4][1])) +
                                   v * (corners[6][1] + u * (corners[7][1] - corners[6][1])));
                const FT pz = (1 - w) * ((1 - v) * (corners[0][2] + u * (corners[1][2] - corners[0][2])) +
                                         v * (corners[2][2] + u * (corners[3][2] - corners[2][2]))) +
                              w * ((1 - v) * (corners[4][2] + u * (corners[5][2] - corners[4][2])) +
                                   v * (corners[6][2] + u * (corners[7][2] - corners[6][2])));

                tg_idx[i] = (uint)points.size();
                points.push_back(Point(px, py, pz));
            }

            // triangulate contours with inner hexagon
            unsigned char tcon_[12];
            for (int i = 0; i < (int)cnt_; i++) {
                if (_not_tunnel != i) {  // contour belongs to tunnel
                    const int cnt_sz = (int)get_cnt_size(i, c_);
                    for (int r = 0; r < cnt_sz; r++) {
                        uint index = -1;
                        double dist = 1000.;
                        uint ci = get_c(i, r, c_);
                        const double u_edge = e_vert(ci, 0);
                        const double v_edge = e_vert(ci, 1);
                        const double w_edge = e_vert(ci, 2);
                        for (int s = 0; s < 6; s++) {
                            const double uval = u_edge - hvt[s][0];
                            const double vval = v_edge - hvt[s][1];
                            const double wval = w_edge - hvt[s][2];
                            double val = uval * uval + vval * vval + wval * wval;
                            if (dist > val) {
                                index = s;
                                dist = val;
                            }
                        }
                        tcon_[ci] = (unsigned char)index;
                    }

                    // correspondence between vertices found
                    // create triangles
                    // needs some functions
                    auto distanceRingIntsModulo = [](const int d1, const int d2) {
                        const int r = (d1 - d2) < 0 ? d2 - d1 : d1 - d2;
                        return (r > 2 ? 6 - r : r);
                    };
                    auto midpointRingIntModulo = [](const int d1, const int d2) {
                        const int dmax = (d1 > d2) ? d1 : d2;
                        const int dmin = (d1 < d2) ? d1 : d2;
                        return ((dmax + 2) % 6 == dmin) ? (dmax + 1) % 6 : (dmax + dmin) / 2;
                    };

                    for (int r = 0; r < cnt_sz; r++) {
                        const uint tid1 = get_c(i, r, c_);
                        const uint tid2 = get_c(i, ((r + 1) % cnt_sz), c_);
                        const uint cid1 = tcon_[tid1];
                        const uint cid2 = tcon_[tid2];
                        // compute index distance
                        const int dst = distanceRingIntsModulo(cid1, cid2);
                        switch (dst) {
                            case 0: {
                                add_triangle(vertices[tid1], vertices[tid2], tg_idx[cid1]);
                            } break;
                            case 1: {
                                // measure diagonals
                                // triangulate along shortest diagonal
                                double u_edge = e_vert(tid1, 0);
                                double v_edge = e_vert(tid1, 1);
                                double w_edge = e_vert(tid1, 2);
                                const double l1 = (u_edge - hvt[cid2][0]) * (u_edge - hvt[cid2][0]) +
                                                  (v_edge - hvt[cid2][1]) * (v_edge - hvt[cid2][1]) +
                                                  (w_edge - hvt[cid2][2]) * (w_edge - hvt[cid2][2]);
                                u_edge = e_vert(tid2, 0);
                                v_edge = e_vert(tid2, 1);
                                w_edge = e_vert(tid2, 2);
                                const double l2 = (u_edge - hvt[cid1][0]) * (u_edge - hvt[cid1][0]) +
                                                  (v_edge - hvt[cid1][1]) * (v_edge - hvt[cid1][1]) +
                                                  (w_edge - hvt[cid1][2]) * (w_edge - hvt[cid1][2]);

                                if (l1 < l2) {
                                    add_triangle(vertices[tid1], vertices[tid2], tg_idx[cid2]);
                                    add_triangle(vertices[tid1], tg_idx[cid2], tg_idx[cid1]);
                                } else {
                                    add_triangle(vertices[tid1], vertices[tid2], tg_idx[cid1]);
                                    add_triangle(vertices[tid2], tg_idx[cid2], tg_idx[cid1]);
                                }
                            } break;
                            case 2: {
                                const int cidm = midpointRingIntModulo(cid1, cid2);

                                add_triangle(vertices[tid1], vertices[tid2], tg_idx[cidm]);
                                add_triangle(vertices[tid1], tg_idx[cidm], tg_idx[cid1]);
                                add_triangle(vertices[tid2], tg_idx[cid2], tg_idx[cidm]);
                            } break;
                        }  // switch
                    }      // for loop over the vertices of the contour
                }          // if (_not_tunnel)
            }              // for loop over contours
            if (cnt_ == 1) {
                // there is a single contour
                // triangulate and close inner hexagon
                // triangle must have the correct orientation
                // use asymptotic_decider() to see if positive vertices
                // are separated, in thic case orientation must be changed
                const bool s_ = (asymptotic_decider(values[0], values[1], values[2], values[3]) <= i0);
                const bool of_ = (wi[1] < wi[0]) ? s_ : !s_;

                if (!of_) {
                    add_triangle(tg_idx[0], tg_idx[2], tg_idx[1]);
                    add_triangle(tg_idx[2], tg_idx[4], tg_idx[3]);
                    add_triangle(tg_idx[0], tg_idx[5], tg_idx[4]);
                    add_triangle(tg_idx[0], tg_idx[4], tg_idx[2]);
                } else {
                    add_triangle(tg_idx[0], tg_idx[1], tg_idx[2]);
                    add_triangle(tg_idx[2], tg_idx[3], tg_idx[4]);
                    add_triangle(tg_idx[0], tg_idx[4], tg_idx[5]);
                    add_triangle(tg_idx[0], tg_idx[2], tg_idx[4]);
                }
            }
        } else {
            // there is no tunnel
            // handle case with no saddle point as simple polygons with 3, 4, 5 or six vertices
            const unsigned char nr_u{(unsigned char)nrQSolFace(0, q_sol)};
            const unsigned char nr_v{(unsigned char)nrQSolFace(1, q_sol)};
            const unsigned char nr_w{(unsigned char)nrQSolFace(2, q_sol)};
            const unsigned char nr_t{(unsigned char)(nr_u + nr_v + nr_w)};
            if (nr_t == nr_u || nr_t == nr_v || nr_t == nr_w) {
                // loop over all contours
                for (int i = 0; i < (int)cnt_; i++) {
                    switch (get_cnt_size(i, c_)) {
                        case 3: {
                            add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 1, c_)],
                                         vertices[get_c(i, 2, c_)]);
                        } break;
                        case 4: {
                            add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 1, c_)],
                                         vertices[get_c(i, 2, c_)]);
                            add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 2, c_)],
                                         vertices[get_c(i, 3, c_)]);
                        } break;
                        case 5: {
                            add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 1, c_)],
                                         vertices[get_c(i, 2, c_)]);
                            add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 2, c_)],
                                         vertices[get_c(i, 3, c_)]);
                            add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 3, c_)],
                                         vertices[get_c(i, 4, c_)]);
                        } break;
                        case 6: {
                            add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 1, c_)],
                                         vertices[get_c(i, 3, c_)]);
                            add_triangle(vertices[get_c(i, 1, c_)], vertices[get_c(i, 2, c_)],
                                         vertices[get_c(i, 3, c_)]);
                            add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 3, c_)],
                                         vertices[get_c(i, 4, c_)]);
                            add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 4, c_)],
                                         vertices[get_c(i, 5, c_)]);
                        } break;
                    }  // switch over size of contour
                }      // loop over contorus
            }          // thre are no saddle points
            else {
                // there are saddle points
                // fc1 = fs(1, 1)*fs(2, 1) + fs(1, 2)*fs(2, 2);
                // fc2 = fs(1, 1)*fs(3, 1) + fs(1, 2)*fs(3, 2);
                // fc3 = fs(2, 1)*fs(3, 2) + fs(2, 2)*fs(3, 1);
                typedef unsigned char uchar;  // TODO

                unsigned char fs[3][2]{{(uchar)(q_sol & 1), (uchar)((q_sol >> 1) & 1)},
                                       {(uchar)((q_sol >> 2) & 1), (uchar)((q_sol >> 3) & 1)},
                                       {(uchar)((q_sol >> 4) & 1), (uchar)((q_sol >> 5) & 1)}};

                const unsigned char fc1 = fs[0][0] * fs[1][0] + fs[0][1] * fs[1][1];
                const unsigned char fc2 = fs[0][0] * fs[2][0] + fs[0][1] * fs[2][1];
                const unsigned char fc3 = fs[1][0] * fs[2][1] + fs[1][1] * fs[2][0];
                const unsigned char c_faces = fc1 + fc2 + fc3;
                double ucoord{};
                double vcoord{};
                double wcoord{};
                switch (c_faces) {
                    case 2: {
                        if (fc1 == 0) {
                            ucoord = fs[0][0] * ui[0] + fs[0][1] * ui[1];
                            vcoord = fs[1][0] * vi[0] + fs[1][1] * vi[1];
                            wcoord = fs[1][0] * wi[1] + fs[1][1] * wi[0];
                        } else if (fc2 == 0) {
                            ucoord = fs[0][0] * ui[0] + fs[0][1] * ui[1];
                            vcoord = fs[0][0] * vi[0] + fs[0][1] * vi[1];
                            wcoord = fs[0][0] * wi[1] + fs[0][1] * wi[0];
                        } else if (fc3 == 0) {
                            ucoord = fs[1][0] * ui[0] + fs[1][1] * ui[1];
                            vcoord = fs[1][0] * vi[0] + fs[1][1] * vi[1];
                            wcoord = fs[1][0] * wi[0] + fs[1][1] * wi[1];
                        }
                    } break;
                    case 3: {
                        ucoord = (fs[0][0] * ui[0] + fs[0][1] * ui[1]) / (fs[0][0] + fs[0][1]);
                        vcoord = (fs[1][0] * vi[0] + fs[1][1] * vi[1]) / (fs[1][0] + fs[1][1]);
                        wcoord = (fs[2][0] * wi[0] + fs[2][1] * wi[1]) / (fs[2][0] + fs[2][1]);
                    } break;
                    case 4: {
                        const unsigned char nr_u = fs[0][0] + fs[0][1];
                        const unsigned char nr_v = fs[1][0] + fs[1][1];
                        const unsigned char nr_w = fs[2][0] + fs[2][1];
                        if (nr_w == 1) {
                            ucoord = fs[2][0] * ui[0] + fs[2][1] * ui[1];
                            vcoord = fs[2][1] * vi[0] + fs[2][0] * vi[1];
                            wcoord = fs[2][0] * wi[0] + fs[2][1] * wi[1];
                        } else if (nr_v == 1) {
                            ucoord = fs[1][0] * ui[0] + fs[1][1] * ui[1];
                            vcoord = fs[1][0] * vi[0] + fs[1][1] * vi[1];
                            wcoord = fs[1][1] * wi[0] + fs[1][0] * wi[1];
                        } else if (nr_u == 1) {
                            ucoord = fs[0][0] * ui[0] + fs[0][1] * ui[1];
                            vcoord = fs[0][0] * vi[0] + fs[0][1] * vi[1];
                            wcoord = fs[0][0] * wi[0] + fs[0][1] * wi[1];
                        }
                    } break;
                }  // switch(c_faces)

                // create inner vertex
                const FT px =
                    (1 - wcoord) * ((1 - vcoord) * (corners[0][0] + ucoord * (corners[1][0] - corners[0][0])) +
                                    vcoord * (corners[2][0] + ucoord * (corners[3][0] - corners[2][0]))) +
                    wcoord * ((1 - vcoord) * (corners[4][0] + ucoord * (corners[5][0] - corners[4][0])) +
                              vcoord * (corners[6][0] + ucoord * (corners[7][0] - corners[6][0])));
                const FT py =
                    (1 - wcoord) * ((1 - vcoord) * (corners[0][1] + ucoord * (corners[1][1] - corners[0][1])) +
                                    vcoord * (corners[2][1] + ucoord * (corners[3][1] - corners[2][1]))) +
                    wcoord * ((1 - vcoord) * (corners[4][1] + ucoord * (corners[5][1] - corners[4][1])) +
                              vcoord * (corners[6][1] + ucoord * (corners[7][1] - corners[6][1])));
                const FT pz =
                    (1 - wcoord) * ((1 - vcoord) * (corners[0][2] + ucoord * (corners[1][2] - corners[0][2])) +
                                    vcoord * (corners[2][2] + ucoord * (corners[3][2] - corners[2][2]))) +
                    wcoord * ((1 - vcoord) * (corners[4][2] + ucoord * (corners[5][2] - corners[4][2])) +
                              vcoord * (corners[6][2] + ucoord * (corners[7][2] - corners[6][2])));

                const uint g_index = (uint)points.size();
                // loop over the contorus
                bool pt_used = false;
                for (int i = 0; i < (int)cnt_; i++) {
                    const unsigned char cnt_sz = (unsigned char)get_cnt_size(i, c_);
                    if (cnt_sz == 3) {
                        add_triangle(vertices[get_c(i, 0, c_)], vertices[get_c(i, 1, c_)], vertices[get_c(i, 2, c_)]);
                    } else {
                        pt_used = true;
                        for (int t = 0; t < cnt_sz; t++) {
                            add_triangle(vertices[get_c(i, t, c_)], vertices[get_c(i, (t + 1) % cnt_sz, c_)], g_index);
                        }
                    }
                }
                if (pt_used) {
                    points.push_back(Point(px, py, pz));
                }
            }  // else - there are saddle points
        }
    }

private:
    const Domain& domain;
    FT isovalue;

    Point_range& points;
    Polygon_range& polygons;

    // compute a unique global index for vertices
    // use as key the unique edge number
    std::map<Edge_descriptor, std::size_t> vertex_map;

    std::mutex mutex;
};

}  // namespace internal
}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_TMC_INTERNAL_TMC_H
