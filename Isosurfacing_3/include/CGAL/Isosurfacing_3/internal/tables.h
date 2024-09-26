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

#ifndef CGAL_ISOSURFACING_3_INTERNAL_TABLES_H
#define CGAL_ISOSURFACING_3_INTERNAL_TABLES_H

#include <CGAL/license/Isosurfacing_3.h>

namespace CGAL {
namespace Isosurfacing {
namespace internal {
namespace Cube_table {
/*
 * Naming convention from "A parallel dual marching cubes approach
 * to quad only surface reconstruction - Grosso & Zint"
 *
 *        ^ y
 *        |
 *       v2------e2------v3
 *       /|             /|
 *     e11|           e10|
 *     /  e3          /  e1
 *   v6------e6------v7  |
 *    |   |          |   |
 *    |  v0------e0--|---v1 --> x
 *    e7 /           e5 /
 *    | e8           | e9
 *    |/             |/
 *   v4------e4------v5
 *   /
 *  < z
 */

constexpr int N_VERTICES = 8;
constexpr int N_EDGES = 12;

// This table iterates around an edge of a voxel in positive direction, starting from the given voxel (0,0,0). The
// iteration is described in coordinates relative to the given voxel. The last number is the local edge index.
constexpr int edge_to_voxel_neighbor[N_EDGES][4][4] =
{
  {{0, 0, 0, 0}, {0, -1, 0, 2}, {0, -1, -1, 6}, {0, 0, -1, 4}},    // e0
  {{0, 0, 0, 1}, {1, 0, 0, 3}, {1, 0, -1, 7}, {0, 0, -1, 5}},      // e1
  {{0, 0, 0, 2}, {0, 0, -1, 6}, {0, 1, -1, 4}, {0, 1, 0, 0}},      // e2
  {{0, 0, 0, 3}, {0, 0, -1, 7}, {-1, 0, -1, 5}, {-1, 0, 0, 1}},    // e3
  {{0, 0, 0, 4}, {0, 0, 1, 0}, {0, -1, 1, 2}, {0, -1, 0, 6}},      // e4
  {{0, 0, 0, 5}, {0, 0, 1, 1}, {1, 0, 1, 3}, {1, 0, 0, 7}},        // e5
  {{0, 0, 0, 6}, {0, 1, 0, 4}, {0, 1, 1, 0}, {0, 0, 1, 2}},        // e6
  {{0, 0, 0, 7}, {-1, 0, 0, 5}, {-1, 0, 1, 1}, {0, 0, 1, 3}},      // e7
  {{0, 0, 0, 8}, {-1, 0, 0, 9}, {-1, -1, 0, 10}, {0, -1, 0, 11}},  // e8
  {{0, 0, 0, 9}, {0, -1, 0, 10}, {1, -1, 0, 11}, {1, 0, 0, 8}},    // e9
  {{0, 0, 0, 10}, {1, 0, 0, 11}, {1, 1, 0, 8}, {0, 1, 0, 9}},      // e10
  {{0, 0, 0, 11}, {0, 1, 0, 8}, {-1, 1, 0, 9}, {-1, 0, 0, 10}}     // e11
};

/* The global edge index consists of the lexicographical index of the v0 vertex of a voxel, and an index that
 * represents the axis. This table maps from the axis index to the local edge index: 0 = x-axis --> 0 1 = y-axis -->
 * 3 2 = z-axis --> 8
 */
constexpr int edge_store_index[3] = {0, 3, 8};

// The local vertex indices of an edge. The indices are sorted by axis direction.
constexpr int edge_to_vertex[N_EDGES][2] =
{
  {0, 1},  // e0
  {1, 3},  // e1
  {2, 3},  // e2
  {0, 2},  // e3
  {4, 5},  // e4
  {5, 7},  // e5
  {6, 7},  // e6
  {4, 6},  // e7
  {0, 4},  // e8
  {1, 5},  // e9
  {3, 7},  // e10
  {2, 6}   // e11
};

// The local vertex coordinates within a voxel.
constexpr int local_vertex_position[N_VERTICES][3] =
{
  {0, 0, 0},  // v0
  {1, 0, 0},  // v1
  {0, 1, 0},  // v2
  {1, 1, 0},  // v3
  {0, 0, 1},  // v4
  {1, 0, 1},  // v5
  {0, 1, 1},  // v6
  {1, 1, 1}   // v7
};

// Edges are uniquely characterized by the two end vertices, which have a unique vertex id
// the end vertices of the edge are computed in the cell by giving the indices (i,j,k).
// These indices are obtained from the cell index by adding 0 or 1 to i, j or k respectively
// Example: edge 0: (i,j,k) - (i+1,j,k)
//          edge 1: (i+1,j,k) - (i+1,j+1,k)
// The first 3 indices are for the first vertex and the second 3 for the second vertex.
// there are 12 edges, assign to each vertex three edges, the global edge numbering
// consist of 3*global_vertex_id + edge_offset.
constexpr int global_edge_id[][4] = {{0, 0, 0, 0}, {1, 0, 0, 1}, {0, 1, 0, 0}, {0, 0, 0, 1},
                                     {0, 0, 1, 0}, {1, 0, 1, 1}, {0, 1, 1, 0}, {0, 0, 1, 1},
                                     {0, 0, 0, 2}, {1, 0, 0, 2}, {1, 1, 0, 2}, {0, 1, 0, 2}};

// probably a list without errors
// indicates which edges has to be intersected
constexpr int intersected_edges[256] =
{
  0,    265,  515,  778,  2060, 2309, 2575, 2822, 1030, 1295, 1541, 1804, 3082, 3331, 3593, 3840, 400,  153,  915,
  666,  2460, 2197, 2975, 2710, 1430, 1183, 1941, 1692, 3482, 3219, 3993, 3728, 560,  825,  51,   314,  2620, 2869,
  2111, 2358, 1590, 1855, 1077, 1340, 3642, 3891, 3129, 3376, 928,  681,  419,  170,  2988, 2725, 2479, 2214, 1958,
  1711, 1445, 1196, 4010, 3747, 3497, 3232, 2240, 2505, 2755, 3018, 204,  453,  719,  966,  3270, 3535, 3781, 4044,
  1226, 1475, 1737, 1984, 2384, 2137, 2899, 2650, 348,  85,   863,  598,  3414, 3167, 3925, 3676, 1370, 1107, 1881,
  1616, 2800, 3065, 2291, 2554, 764,  1013, 255,  502,  3830, 4095, 3317, 3580, 1786, 2035, 1273, 1520, 2912, 2665,
  2403, 2154, 876,  613,  367,  102,  3942, 3695, 3429, 3180, 1898, 1635, 1385, 1120, 1120, 1385, 1635, 1898, 3180,
  3429, 3695, 3942, 102,  367,  613,  876,  2154, 2403, 2665, 2912, 1520, 1273, 2035, 1786, 3580, 3317, 4095, 3830,
  502,  255,  1013, 764,  2554, 2291, 3065, 2800, 1616, 1881, 1107, 1370, 3676, 3925, 3167, 3414, 598,  863,  85,
  348,  2650, 2899, 2137, 2384, 1984, 1737, 1475, 1226, 4044, 3781, 3535, 3270, 966,  719,  453,  204,  3018, 2755,
  2505, 2240, 3232, 3497, 3747, 4010, 1196, 1445, 1711, 1958, 2214, 2479, 2725, 2988, 170,  419,  681,  928,  3376,
  3129, 3891, 3642, 1340, 1077, 1855, 1590, 2358, 2111, 2869, 2620, 314,  51,   825,  560,  3728, 3993, 3219, 3482,
  1692, 1941, 1183, 1430, 2710, 2975, 2197, 2460, 666,  915,  153,  400,  3840, 3593, 3331, 3082, 1804, 1541, 1295,
  1030, 2822, 2575, 2309, 2060, 778,  515,  265,  0
};

// list of triangles for Marching Cubes case t, position at t*16 + tri
constexpr int triangle_cases[4096] =
{
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 0 <-> mc: 0, class rep: 0
  0,  3,  8,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 1 <-> mc: 1, class rep: 1
  0,  9,  1,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 2 <-> mc: 2, class rep: 1
  1,  3,  8,  9,  1,  8,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 3 <-> mc: 3, class rep: 3
  3,  2,  11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 4 <-> mc: 8, class rep: 1
  0,  2,  11, 8,  0,  11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 5 <-> mc: 9, class rep: 3
  1,  0,  9,  2,  11, 3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 6 <-> mc: 10, class rep: 6
  1,  2,  11, 1,  11, 9,  9,  11, 8,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 7 <-> mc: 11, class rep: 7
  1,  10, 2,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 8 <-> mc: 4, class rep: 1
  0,  3,  8,  1,  10, 2,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 9 <-> mc: 5, class rep: 6
  9,  10, 2,  0,  9,  2,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 10 <-> mc: 6, class rep: 3
  2,  3,  8,  2,  8,  10, 10, 8,  9,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 11 <-> mc: 7, class rep: 7
  3,  1,  10, 11, 3,  10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 12 <-> mc: 12, class rep: 3
  0,  1,  10, 0,  10, 8,  8,  10, 11, -1, -1, -1, -1, -1, -1, -1,  // quitte: 13 <-> mc: 13, class rep: 7
  3,  0,  9,  3,  9,  11, 11, 9,  10, -1, -1, -1, -1, -1, -1, -1,  // quitte: 14 <-> mc: 14, class rep: 7
  9,  10, 8,  10, 11, 8,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 15 <-> mc: 15, class rep: 15
  4,  8,  7,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 16 <-> mc: 16, class rep: 1
  4,  0,  3,  7,  4,  3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 17 <-> mc: 17, class rep: 3
  0,  9,  1,  8,  7,  4,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 18 <-> mc: 18, class rep: 6
  4,  9,  1,  4,  1,  7,  7,  1,  3,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 19 <-> mc: 19, class rep: 7
  8,  7,  4,  3,  2,  11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 20 <-> mc: 24, class rep: 6
  11, 7,  4,  11, 4,  2,  2,  4,  0,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 21 <-> mc: 25, class rep: 7
  9,  1,  0,  8,  7,  4,  2,  11, 3,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 22 <-> mc: 26, class rep: 22
  4,  11, 7,  9,  11, 4,  9,  2,  11, 9,  1,  2,  -1, -1, -1, -1,  // quitte: 23 <-> mc: 27, class rep: 23
  1,  10, 2,  8,  7,  4,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 24 <-> mc: 20, class rep: 24
  3,  7,  4,  3,  4,  0,  1,  10, 2,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 25 <-> mc: 21, class rep: 25
  9,  10, 2,  9,  2,  0,  8,  7,  4,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 26 <-> mc: 22, class rep: 25
  2,  9,  10, 2,  7,  9,  2,  3,  7,  7,  4,  9,  -1, -1, -1, -1,  // quitte: 27 <-> mc: 23, class rep: 27
  3,  1,  10, 3,  10, 11, 7,  4,  8,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 28 <-> mc: 28, class rep: 25
  1,  10, 11, 1,  11, 4,  1,  4,  0,  7,  4,  11, -1, -1, -1, -1,  // quitte: 29 <-> mc: 29, class rep: 29
  4,  8,  7,  9,  11, 0,  9,  10, 11, 11, 3,  0,  -1, -1, -1, -1,  // quitte: 30 <-> mc: 30, class rep: 30
  4,  11, 7,  4,  9,  11, 9,  10, 11, -1, -1, -1, -1, -1, -1, -1,  // quitte: 31 <-> mc: 31, class rep: 7
  9,  4,  5,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 32 <-> mc: 32, class rep: 1
  9,  4,  5,  0,  3,  8,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 33 <-> mc: 33, class rep: 6
  0,  4,  5,  1,  0,  5,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 34 <-> mc: 34, class rep: 3
  8,  4,  5,  8,  5,  3,  3,  5,  1,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 35 <-> mc: 35, class rep: 7
  9,  4,  5,  2,  11, 3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 36 <-> mc: 40, class rep: 24
  0,  2,  11, 0,  11, 8,  4,  5,  9,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 37 <-> mc: 41, class rep: 25
  0,  4,  5,  0,  5,  1,  2,  11, 3,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 38 <-> mc: 42, class rep: 25
  2,  5,  1,  2,  8,  5,  2,  11, 8,  4,  5,  8,  -1, -1, -1, -1,  // quitte: 39 <-> mc: 43, class rep: 29
  1,  10, 2,  9,  4,  5,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 40 <-> mc: 36, class rep: 6
  3,  8,  0,  1,  10, 2,  4,  5,  9,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 41 <-> mc: 37, class rep: 22
  5,  10, 2,  5,  2,  4,  4,  2,  0,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 42 <-> mc: 38, class rep: 7
  2,  5,  10, 3,  5,  2,  3,  4,  5,  3,  8,  4,  -1, -1, -1, -1,  // quitte: 43 <-> mc: 39, class rep: 23
  10, 11, 3,  10, 3,  1,  9,  4,  5,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 44 <-> mc: 44, class rep: 25
  4,  5,  9,  0,  1,  8,  8,  1,  10, 8,  10, 11, -1, -1, -1, -1,  // quitte: 45 <-> mc: 45, class rep: 30
  5,  0,  4,  5,  11, 0,  5,  10, 11, 11, 3,  0,  -1, -1, -1, -1,  // quitte: 46 <-> mc: 46, class rep: 27
  5,  8,  4,  5,  10, 8,  10, 11, 8,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 47 <-> mc: 47, class rep: 7
  9,  8,  7,  5,  9,  7,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 48 <-> mc: 48, class rep: 3
  9,  0,  3,  9,  3,  5,  5,  3,  7,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 49 <-> mc: 49, class rep: 7
  0,  8,  7,  0,  7,  1,  1,  7,  5,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 50 <-> mc: 50, class rep: 7
  1,  3,  5,  3,  7,  5,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 51 <-> mc: 51, class rep: 15
  7,  5,  9,  7,  9,  8,  3,  2,  11, -1, -1, -1, -1, -1, -1, -1,  // quitte: 52 <-> mc: 56, class rep: 25
  9,  7,  5,  9,  2,  7,  9,  0,  2,  2,  11, 7,  -1, -1, -1, -1,  // quitte: 53 <-> mc: 57, class rep: 27
  2,  11, 3,  0,  8,  1,  1,  8,  7,  1,  7,  5,  -1, -1, -1, -1,  // quitte: 54 <-> mc: 58, class rep: 30
  11, 1,  2,  11, 7,  1,  7,  5,  1,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 55 <-> mc: 59, class rep: 7
  9,  8,  7,  9,  7,  5,  10, 2,  1,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 56 <-> mc: 52, class rep: 25
  10, 2,  1,  9,  0,  5,  5,  0,  3,  5,  3,  7,  -1, -1, -1, -1,  // quitte: 57 <-> mc: 53, class rep: 30
  8,  2,  0,  8,  5,  2,  8,  7,  5,  10, 2,  5,  -1, -1, -1, -1,  // quitte: 58 <-> mc: 54, class rep: 29
  2,  5,  10, 2,  3,  5,  3,  7,  5,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 59 <-> mc: 55, class rep: 7
  9,  8,  5,  8,  7,  5,  10, 3,  1,  10, 11, 3,  -1, -1, -1, -1,  // quitte: 60 <-> mc: 60, class rep: 60
  5,  0,  7,  5,  9,  0,  7,  0,  11, 1,  10, 0,  11, 0,  10, -1,  // quitte: 61 <-> mc: 61, class rep: 25
  11, 0,  10, 11, 3,  0,  10, 0,  5,  8,  7,  0,  5,  0,  7,  -1,  // quitte: 62 <-> mc: 62, class rep: 25
  11, 5,  10, 7,  5,  11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 63 <-> mc: 63, class rep: 3
  7,  11, 6,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 64 <-> mc: 128, class rep: 1
  3,  8,  0,  11, 6,  7,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 65 <-> mc: 129, class rep: 6
  0,  9,  1,  11, 6,  7,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 66 <-> mc: 130, class rep: 24
  8,  9,  1,  8,  1,  3,  11, 6,  7,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 67 <-> mc: 131, class rep: 25
  7,  3,  2,  6,  7,  2,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 68 <-> mc: 136, class rep: 3
  7,  8,  0,  7,  0,  6,  6,  0,  2,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 69 <-> mc: 137, class rep: 7
  2,  6,  7,  2,  7,  3,  0,  9,  1,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 70 <-> mc: 138, class rep: 25
  1,  2,  6,  1,  6,  8,  1,  8,  9,  8,  6,  7,  -1, -1, -1, -1,  // quitte: 71 <-> mc: 139, class rep: 27
  10, 2,  1,  6,  7,  11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 72 <-> mc: 132, class rep: 6
  1,  10, 2,  3,  8,  0,  6,  7,  11, -1, -1, -1, -1, -1, -1, -1,  // quitte: 73 <-> mc: 133, class rep: 22
  2,  0,  9,  2,  9,  10, 6,  7,  11, -1, -1, -1, -1, -1, -1, -1,  // quitte: 74 <-> mc: 134, class rep: 25
  6,  7,  11, 2,  3,  10, 10, 3,  8,  10, 8,  9,  -1, -1, -1, -1,  // quitte: 75 <-> mc: 135, class rep: 30
  10, 6,  7,  10, 7,  1,  1,  7,  3,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 76 <-> mc: 140, class rep: 7
  10, 6,  7,  1,  10, 7,  1,  7,  8,  1,  8,  0,  -1, -1, -1, -1,  // quitte: 77 <-> mc: 141, class rep: 23
  0,  7,  3,  0,  10, 7,  0,  9,  10, 6,  7,  10, -1, -1, -1, -1,  // quitte: 78 <-> mc: 142, class rep: 29
  7,  10, 6,  7,  8,  10, 8,  9,  10, -1, -1, -1, -1, -1, -1, -1,  // quitte: 79 <-> mc: 143, class rep: 7
  6,  4,  8,  11, 6,  8,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 80 <-> mc: 144, class rep: 3
  3,  11, 6,  3,  6,  0,  0,  6,  4,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 81 <-> mc: 145, class rep: 7
  8,  11, 6,  8,  6,  4,  9,  1,  0,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 82 <-> mc: 146, class rep: 25
  9,  6,  4,  9,  3,  6,  9,  1,  3,  11, 6,  3,  -1, -1, -1, -1,  // quitte: 83 <-> mc: 147, class rep: 29
  8,  3,  2,  8,  2,  4,  4,  2,  6,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 84 <-> mc: 152, class rep: 7
  0,  2,  4,  4,  2,  6,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 85 <-> mc: 153, class rep: 15
  1,  0,  9,  2,  4,  3,  2,  6,  4,  4,  8,  3,  -1, -1, -1, -1,  // quitte: 86 <-> mc: 154, class rep: 30
  1,  4,  9,  1,  2,  4,  2,  6,  4,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 87 <-> mc: 155, class rep: 7
  6,  4,  8,  6,  8,  11, 2,  1,  10, -1, -1, -1, -1, -1, -1, -1,  // quitte: 88 <-> mc: 148, class rep: 25
  1,  10, 2,  3,  11, 0,  0,  11, 6,  0,  6,  4,  -1, -1, -1, -1,  // quitte: 89 <-> mc: 149, class rep: 30
  4,  8,  11, 4,  11, 6,  0,  9,  2,  2,  9,  10, -1, -1, -1, -1,  // quitte: 90 <-> mc: 150, class rep: 60
  10, 3,  9,  10, 2,  3,  9,  3,  4,  11, 6,  3,  4,  3,  6,  -1,  // quitte: 91 <-> mc: 151, class rep: 25
  8,  3,  1,  8,  1,  6,  8,  6,  4,  6,  1,  10, -1, -1, -1, -1,  // quitte: 92 <-> mc: 156, class rep: 27
  10, 0,  1,  10, 6,  0,  6,  4,  0,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 93 <-> mc: 157, class rep: 7
  4,  3,  6,  4,  8,  3,  6,  3,  10, 0,  9,  3,  10, 3,  9,  -1,  // quitte: 94 <-> mc: 158, class rep: 25
  10, 4,  9,  6,  4,  10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 95 <-> mc: 159, class rep: 3
  4,  5,  9,  7,  11, 6,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 96 <-> mc: 160, class rep: 6
  0,  3,  8,  4,  5,  9,  11, 6,  7,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 97 <-> mc: 161, class rep: 22
  5,  1,  0,  5,  0,  4,  7,  11, 6,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 98 <-> mc: 162, class rep: 25
  11, 6,  7,  8,  4,  3,  3,  4,  5,  3,  5,  1,  -1, -1, -1, -1,  // quitte: 99 <-> mc: 163, class rep: 30
  7,  3,  2,  7,  2,  6,  5,  9,  4,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 100 <-> mc: 168, class rep: 25
  9,  4,  5,  0,  6,  8,  0,  2,  6,  6,  7,  8,  -1, -1, -1, -1,  // quitte: 101 <-> mc: 169, class rep: 30
  3,  2,  6,  3,  6,  7,  1,  0,  5,  5,  0,  4,  -1, -1, -1, -1,  // quitte: 102 <-> mc: 170, class rep: 60
  6,  8,  2,  6,  7,  8,  2,  8,  1,  4,  5,  8,  1,  8,  5,  -1,  // quitte: 103 <-> mc: 171, class rep: 25
  9,  4,  5,  10, 2,  1,  7,  11, 6,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 104 <-> mc: 164, class rep: 22
  6,  7,  11, 1,  10, 2,  0,  3,  8,  4,  5,  9,  -1, -1, -1, -1,  // quitte: 105 <-> mc: 165, class rep: 105
  7,  11, 6,  5,  10, 4,  4,  10, 2,  4,  2,  0,  -1, -1, -1, -1,  // quitte: 106 <-> mc: 166, class rep: 30
  3,  8,  4,  3,  4,  5,  3,  5,  2,  10, 2,  5,  11, 6,  7,  -1,  // quitte: 107 <-> mc: 167, class rep: 22
  9,  4,  5,  10, 6,  1,  1,  6,  7,  1,  7,  3,  -1, -1, -1, -1,  // quitte: 108 <-> mc: 172, class rep: 30
  1,  10, 6,  1,  6,  7,  1,  7,  0,  8,  0,  7,  9,  4,  5,  -1,  // quitte: 109 <-> mc: 173, class rep: 22
  4,  10, 0,  4,  5,  10, 0,  10, 3,  6,  7,  10, 3,  10, 7,  -1,  // quitte: 110 <-> mc: 174, class rep: 25
  7,  10, 6,  7,  8,  10, 5,  10, 4,  4,  10, 8,  -1, -1, -1, -1,  // quitte: 111 <-> mc: 175, class rep: 6
  6,  5,  9,  6,  9,  11, 11, 9,  8,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 112 <-> mc: 176, class rep: 7
  3,  11, 6,  0,  3,  6,  0,  6,  5,  0,  5,  9,  -1, -1, -1, -1,  // quitte: 113 <-> mc: 177, class rep: 23
  0,  8,  11, 0,  11, 5,  0,  5,  1,  5,  11, 6,  -1, -1, -1, -1,  // quitte: 114 <-> mc: 178, class rep: 27
  6,  3,  11, 6,  5,  3,  5,  1,  3,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 115 <-> mc: 179, class rep: 7
  5,  9,  8,  5,  8,  2,  5,  2,  6,  3,  2,  8,  -1, -1, -1, -1,  // quitte: 116 <-> mc: 184, class rep: 29
  9,  6,  5,  9,  0,  6,  0,  2,  6,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 117 <-> mc: 185, class rep: 7
  1,  8,  5,  1,  0,  8,  5,  8,  6,  3,  2,  8,  6,  8,  2,  -1,  // quitte: 118 <-> mc: 186, class rep: 25
  1,  6,  5,  2,  6,  1,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 119 <-> mc: 187, class rep: 3
  1,  10, 2,  9,  11, 5,  9,  8,  11, 11, 6,  5,  -1, -1, -1, -1,  // quitte: 120 <-> mc: 180, class rep: 30
  0,  3,  11, 0,  11, 6,  0,  6,  9,  5,  9,  6,  1,  10, 2,  -1,  // quitte: 121 <-> mc: 181, class rep: 22
  11, 5,  8,  11, 6,  5,  8,  5,  0,  10, 2,  5,  0,  5,  2,  -1,  // quitte: 122 <-> mc: 182, class rep: 25
  6,  3,  11, 6,  5,  3,  2,  3,  10, 10, 3,  5,  -1, -1, -1, -1,  // quitte: 123 <-> mc: 183, class rep: 6
  1,  6,  3,  1,  10, 6,  3,  6,  8,  5,  9,  6,  8,  6,  9,  -1,  // quitte: 124 <-> mc: 188, class rep: 25
  10, 0,  1,  10, 6,  0,  9,  0,  5,  5,  0,  6,  -1, -1, -1, -1,  // quitte: 125 <-> mc: 189, class rep: 6
  0,  8,  3,  5,  10, 6,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 126 <-> mc: 190, class rep: 24
  10, 6,  5,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 127 <-> mc: 191, class rep: 1
  10, 5,  6,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 128 <-> mc: 64, class rep: 1
  0,  3,  8,  5,  6,  10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 129 <-> mc: 65, class rep: 24
  9,  1,  0,  5,  6,  10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 130 <-> mc: 66, class rep: 6
  1,  3,  8,  1,  8,  9,  5,  6,  10, -1, -1, -1, -1, -1, -1, -1,  // quitte: 131 <-> mc: 67, class rep: 25
  2,  11, 3,  10, 5,  6,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 132 <-> mc: 72, class rep: 6
  11, 8,  0,  11, 0,  2,  10, 5,  6,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 133 <-> mc: 73, class rep: 25
  0,  9,  1,  2,  11, 3,  5,  6,  10, -1, -1, -1, -1, -1, -1, -1,  // quitte: 134 <-> mc: 74, class rep: 22
  5,  6,  10, 1,  2,  9,  9,  2,  11, 9,  11, 8,  -1, -1, -1, -1,  // quitte: 135 <-> mc: 75, class rep: 30
  1,  5,  6,  2,  1,  6,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 136 <-> mc: 68, class rep: 3
  1,  5,  6,  1,  6,  2,  3,  8,  0,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 137 <-> mc: 69, class rep: 25
  9,  5,  6,  9,  6,  0,  0,  6,  2,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 138 <-> mc: 70, class rep: 7
  5,  8,  9,  5,  2,  8,  5,  6,  2,  3,  8,  2,  -1, -1, -1, -1,  // quitte: 139 <-> mc: 71, class rep: 29
  6,  11, 3,  6,  3,  5,  5,  3,  1,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 140 <-> mc: 76, class rep: 7
  0,  11, 8,  0,  5,  11, 0,  1,  5,  5,  6,  11, -1, -1, -1, -1,  // quitte: 141 <-> mc: 77, class rep: 27
  3,  6,  11, 0,  6,  3,  0,  5,  6,  0,  9,  5,  -1, -1, -1, -1,  // quitte: 142 <-> mc: 78, class rep: 23
  6,  9,  5,  6,  11, 9,  11, 8,  9,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 143 <-> mc: 79, class rep: 7
  5,  6,  10, 4,  8,  7,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 144 <-> mc: 80, class rep: 6
  4,  0,  3,  4,  3,  7,  6,  10, 5,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 145 <-> mc: 81, class rep: 25
  1,  0,  9,  5,  6,  10, 8,  7,  4,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 146 <-> mc: 82, class rep: 22
  10, 5,  6,  1,  7,  9,  1,  3,  7,  7,  4,  9,  -1, -1, -1, -1,  // quitte: 147 <-> mc: 83, class rep: 30
  3,  2,  11, 7,  4,  8,  10, 5,  6,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 148 <-> mc: 88, class rep: 22
  5,  6,  10, 4,  2,  7,  4,  0,  2,  2,  11, 7,  -1, -1, -1, -1,  // quitte: 149 <-> mc: 89, class rep: 30
  0,  9,  1,  4,  8,  7,  2,  11, 3,  5,  6,  10, -1, -1, -1, -1,  // quitte: 150 <-> mc: 90, class rep: 105
  9,  1,  2,  9,  2,  11, 9,  11, 4,  7,  4,  11, 5,  6,  10, -1,  // quitte: 151 <-> mc: 91, class rep: 22
  6,  2,  1,  6,  1,  5,  4,  8,  7,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 152 <-> mc: 84, class rep: 25
  1,  5,  2,  5,  6,  2,  3,  4,  0,  3,  7,  4,  -1, -1, -1, -1,  // quitte: 153 <-> mc: 85, class rep: 60
  8,  7,  4,  9,  5,  0,  0,  5,  6,  0,  6,  2,  -1, -1, -1, -1,  // quitte: 154 <-> mc: 86, class rep: 30
  7,  9,  3,  7,  4,  9,  3,  9,  2,  5,  6,  9,  2,  9,  6,  -1,  // quitte: 155 <-> mc: 87, class rep: 25
  8,  7,  4,  3,  5,  11, 3,  1,  5,  5,  6,  11, -1, -1, -1, -1,  // quitte: 156 <-> mc: 92, class rep: 30
  5,  11, 1,  5,  6,  11, 1,  11, 0,  7,  4,  11, 0,  11, 4,  -1,  // quitte: 157 <-> mc: 93, class rep: 25
  0,  9,  5,  0,  5,  6,  0,  6,  3,  11, 3,  6,  8,  7,  4,  -1,  // quitte: 158 <-> mc: 94, class rep: 22
  6,  9,  5,  6,  11, 9,  4,  9,  7,  7,  9,  11, -1, -1, -1, -1,  // quitte: 159 <-> mc: 95, class rep: 6
  10, 9,  4,  6,  10, 4,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 160 <-> mc: 96, class rep: 3
  4,  6,  10, 4,  10, 9,  0,  3,  8,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 161 <-> mc: 97, class rep: 25
  10, 1,  0,  10, 0,  6,  6,  0,  4,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 162 <-> mc: 98, class rep: 7
  8,  1,  3,  8,  6,  1,  8,  4,  6,  6,  10, 1,  -1, -1, -1, -1,  // quitte: 163 <-> mc: 99, class rep: 27
  10, 9,  4,  10, 4,  6,  11, 3,  2,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 164 <-> mc: 104, class rep: 25
  0,  2,  8,  2,  11, 8,  4,  10, 9,  4,  6,  10, -1, -1, -1, -1,  // quitte: 165 <-> mc: 105, class rep: 60
  3,  2,  11, 0,  6,  1,  0,  4,  6,  6,  10, 1,  -1, -1, -1, -1,  // quitte: 166 <-> mc: 106, class rep: 30
  6,  1,  4,  6,  10, 1,  4,  1,  8,  2,  11, 1,  8,  1,  11, -1,  // quitte: 167 <-> mc: 107, class rep: 25
  1,  9,  4,  1,  4,  2,  2,  4,  6,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 168 <-> mc: 100, class rep: 7
  3,  8,  0,  1,  9,  2,  2,  9,  4,  2,  4,  6,  -1, -1, -1, -1,  // quitte: 169 <-> mc: 101, class rep: 30
  0,  4,  2,  4,  6,  2,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 170 <-> mc: 102, class rep: 15
  8,  2,  3,  8,  4,  2,  4,  6,  2,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 171 <-> mc: 103, class rep: 7
  9,  4,  6,  9,  6,  3,  9,  3,  1,  11, 3,  6,  -1, -1, -1, -1,  // quitte: 172 <-> mc: 108, class rep: 29
  8,  1,  11, 8,  0,  1,  11, 1,  6,  9,  4,  1,  6,  1,  4,  -1,  // quitte: 173 <-> mc: 109, class rep: 25
  3,  6,  11, 3,  0,  6,  0,  4,  6,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 174 <-> mc: 110, class rep: 7
  6,  8,  4,  11, 8,  6,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 175 <-> mc: 111, class rep: 3
  7,  6,  10, 7,  10, 8,  8,  10, 9,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 176 <-> mc: 112, class rep: 7
  0,  3,  7,  0,  7,  10, 0,  10, 9,  6,  10, 7,  -1, -1, -1, -1,  // quitte: 177 <-> mc: 113, class rep: 29
  10, 7,  6,  1,  7,  10, 1,  8,  7,  1,  0,  8,  -1, -1, -1, -1,  // quitte: 178 <-> mc: 114, class rep: 23
  10, 7,  6,  10, 1,  7,  1,  3,  7,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 179 <-> mc: 115, class rep: 7
  2,  11, 3,  10, 8,  6,  10, 9,  8,  8,  7,  6,  -1, -1, -1, -1,  // quitte: 180 <-> mc: 120, class rep: 30
  2,  7,  0,  2,  11, 7,  0,  7,  9,  6,  10, 7,  9,  7,  10, -1,  // quitte: 181 <-> mc: 121, class rep: 25
  1,  0,  8,  1,  8,  7,  1,  7,  10, 6,  10, 7,  2,  11, 3,  -1,  // quitte: 182 <-> mc: 122, class rep: 22
  11, 1,  2,  11, 7,  1,  10, 1,  6,  6,  1,  7,  -1, -1, -1, -1,  // quitte: 183 <-> mc: 123, class rep: 6
  1,  6,  2,  1,  8,  6,  1,  9,  8,  8,  7,  6,  -1, -1, -1, -1,  // quitte: 184 <-> mc: 116, class rep: 27
  2,  9,  6,  2,  1,  9,  6,  9,  7,  0,  3,  9,  7,  9,  3,  -1,  // quitte: 185 <-> mc: 117, class rep: 25
  7,  0,  8,  7,  6,  0,  6,  2,  0,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 186 <-> mc: 118, class rep: 7
  7,  2,  3,  6,  2,  7,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 187 <-> mc: 119, class rep: 3
  8,  6,  9,  8,  7,  6,  9,  6,  1,  11, 3,  6,  1,  6,  3,  -1,  // quitte: 188 <-> mc: 124, class rep: 25
  0,  1,  9,  11, 7,  6,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 189 <-> mc: 125, class rep: 24
  7,  0,  8,  7,  6,  0,  3,  0,  11, 11, 0,  6,  -1, -1, -1, -1,  // quitte: 190 <-> mc: 126, class rep: 6
  7,  6,  11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 191 <-> mc: 127, class rep: 1
  11, 10, 5,  7,  11, 5,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 192 <-> mc: 192, class rep: 3
  11, 10, 5,  11, 5,  7,  8,  0,  3,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 193 <-> mc: 193, class rep: 25
  5,  7,  11, 5,  11, 10, 1,  0,  9,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 194 <-> mc: 194, class rep: 25
  10, 5,  7,  10, 7,  11, 9,  1,  8,  8,  1,  3,  -1, -1, -1, -1,  // quitte: 195 <-> mc: 195, class rep: 60
  2,  10, 5,  2,  5,  3,  3,  5,  7,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 196 <-> mc: 200, class rep: 7
  8,  0,  2,  8,  2,  5,  8,  5,  7,  10, 5,  2,  -1, -1, -1, -1,  // quitte: 197 <-> mc: 201, class rep: 29
  9,  1,  0,  5,  3,  10, 5,  7,  3,  3,  2,  10, -1, -1, -1, -1,  // quitte: 198 <-> mc: 202, class rep: 30
  9,  2,  8,  9,  1,  2,  8,  2,  7,  10, 5,  2,  7,  2,  5,  -1,  // quitte: 199 <-> mc: 203, class rep: 25
  11, 2,  1,  11, 1,  7,  7,  1,  5,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 200 <-> mc: 196, class rep: 7
  0,  3,  8,  1,  7,  2,  1,  5,  7,  7,  11, 2,  -1, -1, -1, -1,  // quitte: 201 <-> mc: 197, class rep: 30
  9,  5,  7,  9,  7,  2,  9,  2,  0,  2,  7,  11, -1, -1, -1, -1,  // quitte: 202 <-> mc: 198, class rep: 27
  7,  2,  5,  7,  11, 2,  5,  2,  9,  3,  8,  2,  9,  2,  8,  -1,  // quitte: 203 <-> mc: 199, class rep: 25
  1,  5,  3,  3,  5,  7,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 204 <-> mc: 204, class rep: 15
  0,  7,  8,  0,  1,  7,  1,  5,  7,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 205 <-> mc: 205, class rep: 7
  9,  3,  0,  9,  5,  3,  5,  7,  3,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 206 <-> mc: 206, class rep: 7
  9,  7,  8,  5,  7,  9,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 207 <-> mc: 207, class rep: 3
  5,  4,  8,  5,  8,  10, 10, 8,  11, -1, -1, -1, -1, -1, -1, -1,  // quitte: 208 <-> mc: 208, class rep: 7
  5,  4,  0,  5,  0,  11, 5,  11, 10, 11, 0,  3,  -1, -1, -1, -1,  // quitte: 209 <-> mc: 209, class rep: 27
  0,  9,  1,  8,  10, 4,  8,  11, 10, 10, 5,  4,  -1, -1, -1, -1,  // quitte: 210 <-> mc: 210, class rep: 30
  10, 4,  11, 10, 5,  4,  11, 4,  3,  9,  1,  4,  3,  4,  1,  -1,  // quitte: 211 <-> mc: 211, class rep: 25
  2,  10, 5,  3,  2,  5,  3,  5,  4,  3,  4,  8,  -1, -1, -1, -1,  // quitte: 212 <-> mc: 216, class rep: 23
  5,  2,  10, 5,  4,  2,  4,  0,  2,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 213 <-> mc: 217, class rep: 7
  3,  2,  10, 3,  10, 5,  3,  5,  8,  4,  8,  5,  0,  9,  1,  -1,  // quitte: 214 <-> mc: 218, class rep: 22
  5,  2,  10, 5,  4,  2,  1,  2,  9,  9,  2,  4,  -1, -1, -1, -1,  // quitte: 215 <-> mc: 219, class rep: 6
  2,  1,  5,  2,  5,  8,  2,  8,  11, 4,  8,  5,  -1, -1, -1, -1,  // quitte: 216 <-> mc: 212, class rep: 29
  0,  11, 4,  0,  3,  11, 4,  11, 5,  2,  1,  11, 5,  11, 1,  -1,  // quitte: 217 <-> mc: 213, class rep: 25
  0,  5,  2,  0,  9,  5,  2,  5,  11, 4,  8,  5,  11, 5,  8,  -1,  // quitte: 218 <-> mc: 214, class rep: 25
  9,  5,  4,  2,  3,  11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 219 <-> mc: 215, class rep: 24
  8,  5,  4,  8,  3,  5,  3,  1,  5,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 220 <-> mc: 220, class rep: 7
  0,  5,  4,  1,  5,  0,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 221 <-> mc: 221, class rep: 3
  8,  5,  4,  8,  3,  5,  9,  5,  0,  0,  5,  3,  -1, -1, -1, -1,  // quitte: 222 <-> mc: 222, class rep: 6
  9,  5,  4,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 223 <-> mc: 223, class rep: 1
  4,  7,  11, 4,  11, 9,  9,  11, 10, -1, -1, -1, -1, -1, -1, -1,  // quitte: 224 <-> mc: 224, class rep: 7
  0,  3,  8,  4,  7,  9,  9,  7,  11, 9,  11, 10, -1, -1, -1, -1,  // quitte: 225 <-> mc: 225, class rep: 30
  1,  11, 10, 1,  4,  11, 1,  0,  4,  7,  11, 4,  -1, -1, -1, -1,  // quitte: 226 <-> mc: 226, class rep: 29
  3,  4,  1,  3,  8,  4,  1,  4,  10, 7,  11, 4,  10, 4,  11, -1,  // quitte: 227 <-> mc: 227, class rep: 25
  2,  10, 9,  2,  9,  7,  2,  7,  3,  7,  9,  4,  -1, -1, -1, -1,  // quitte: 228 <-> mc: 232, class rep: 27
  9,  7,  10, 9,  4,  7,  10, 7,  2,  8,  0,  7,  2,  7,  0,  -1,  // quitte: 229 <-> mc: 233, class rep: 25
  3,  10, 7,  3,  2,  10, 7,  10, 4,  1,  0,  10, 4,  10, 0,  -1,  // quitte: 230 <-> mc: 234, class rep: 25
  1,  2,  10, 8,  4,  7,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 231 <-> mc: 235, class rep: 24
  4,  7,  11, 9,  4,  11, 9,  11, 2,  9,  2,  1,  -1, -1, -1, -1,  // quitte: 232 <-> mc: 228, class rep: 23
  9,  4,  7,  9,  7,  11, 9,  11, 1,  2,  1,  11, 0,  3,  8,  -1,  // quitte: 233 <-> mc: 229, class rep: 22
  11, 4,  7,  11, 2,  4,  2,  0,  4,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 234 <-> mc: 230, class rep: 7
  11, 4,  7,  11, 2,  4,  8,  4,  3,  3,  4,  2,  -1, -1, -1, -1,  // quitte: 235 <-> mc: 231, class rep: 6
  4,  1,  9,  4,  7,  1,  7,  3,  1,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 236 <-> mc: 236, class rep: 7
  4,  1,  9,  4,  7,  1,  0,  1,  8,  8,  1,  7,  -1, -1, -1, -1,  // quitte: 237 <-> mc: 237, class rep: 6
  4,  3,  0,  7,  3,  4,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 238 <-> mc: 238, class rep: 3
  4,  7,  8,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 239 <-> mc: 239, class rep: 1
  9,  8,  10, 10, 8,  11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 240 <-> mc: 240, class rep: 15
  3,  9,  0,  3,  11, 9,  11, 10, 9,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 241 <-> mc: 241, class rep: 7
  0,  10, 1,  0,  8,  10, 8,  11, 10, -1, -1, -1, -1, -1, -1, -1,  // quitte: 242 <-> mc: 242, class rep: 7
  3,  10, 1,  11, 10, 3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 243 <-> mc: 243, class rep: 3
  2,  8,  3,  2,  10, 8,  10, 9,  8,  -1, -1, -1, -1, -1, -1, -1,  // quitte: 244 <-> mc: 248, class rep: 7
  9,  2,  10, 0,  2,  9,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 245 <-> mc: 249, class rep: 3
  2,  8,  3,  2,  10, 8,  0,  8,  1,  1,  8,  10, -1, -1, -1, -1,  // quitte: 246 <-> mc: 250, class rep: 6
  1,  2,  10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 247 <-> mc: 251, class rep: 1
  1,  11, 2,  1,  9,  11, 9,  8,  11, -1, -1, -1, -1, -1, -1, -1,  // quitte: 248 <-> mc: 244, class rep: 7
  3,  9,  0,  3,  11, 9,  1,  9,  2,  2,  9,  11, -1, -1, -1, -1,  // quitte: 249 <-> mc: 245, class rep: 6
  0,  11, 2,  8,  11, 0,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 250 <-> mc: 246, class rep: 3
  3,  11, 2,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 251 <-> mc: 247, class rep: 1
  1,  8,  3,  9,  8,  1,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 252 <-> mc: 252, class rep: 3
  0,  1,  9,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 253 <-> mc: 253, class rep: 1
  0,  8,  3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  // quitte: 254 <-> mc: 254, class rep: 1
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1   // quitte: 255 <-> mc: 255, class rep: 0
};

// list of ambiguous cases
constexpr int t_ambig[256] =
{
  0,    // quitte: 0 <-> mc: 0, class representative: 0
  1,    // quitte: 1 <-> mc: 1, class representative: 1
  2,    // quitte: 2 <-> mc: 2, class representative: 1
  3,    // quitte: 3 <-> mc: 3, class representative: 3
  4,    // quitte: 4 <-> mc: 8, class representative: 1
  5,    // quitte: 5 <-> mc: 9, class representative: 3
  105,  // quitte: 6 <-> mc: 10, class representative: 6
  7,    // quitte: 7 <-> mc: 11, class representative: 7
  8,    // quitte: 8 <-> mc: 4, class representative: 1
  105,  // quitte: 9 <-> mc: 5, class representative: 6
  10,   // quitte: 10 <-> mc: 6, class representative: 3
  11,   // quitte: 11 <-> mc: 7, class representative: 7
  12,   // quitte: 12 <-> mc: 12, class representative: 3
  13,   // quitte: 13 <-> mc: 13, class representative: 7
  14,   // quitte: 14 <-> mc: 14, class representative: 7
  15,   // quitte: 15 <-> mc: 15, class representative: 15
  16,   // quitte: 16 <-> mc: 16, class representative: 1
  17,   // quitte: 17 <-> mc: 17, class representative: 3
  105,  // quitte: 18 <-> mc: 18, class representative: 6
  19,   // quitte: 19 <-> mc: 19, class representative: 7
  105,  // quitte: 20 <-> mc: 24, class representative: 6
  21,   // quitte: 21 <-> mc: 25, class representative: 7
  105,  // quitte: 22 <-> mc: 26, class representative: 22
  23,   // quitte: 23 <-> mc: 27, class representative: 23
  105,  // quitte: 24 <-> mc: 20, class representative: 24
  105,  // quitte: 25 <-> mc: 21, class representative: 25
  105,  // quitte: 26 <-> mc: 22, class representative: 25
  27,   // quitte: 27 <-> mc: 23, class representative: 27
  105,  // quitte: 28 <-> mc: 28, class representative: 25
  29,   // quitte: 29 <-> mc: 29, class representative: 29
  105,  // quitte: 30 <-> mc: 30, class representative: 30
  31,   // quitte: 31 <-> mc: 31, class representative: 7
  32,   // quitte: 32 <-> mc: 32, class representative: 1
  105,  // quitte: 33 <-> mc: 33, class representative: 6
  34,   // quitte: 34 <-> mc: 34, class representative: 3
  35,   // quitte: 35 <-> mc: 35, class representative: 7
  105,  // quitte: 36 <-> mc: 40, class representative: 24
  105,  // quitte: 37 <-> mc: 41, class representative: 25
  105,  // quitte: 38 <-> mc: 42, class representative: 25
  39,   // quitte: 39 <-> mc: 43, class representative: 29
  105,  // quitte: 40 <-> mc: 36, class representative: 6
  105,  // quitte: 41 <-> mc: 37, class representative: 22
  42,   // quitte: 42 <-> mc: 38, class representative: 7
  43,   // quitte: 43 <-> mc: 39, class representative: 23
  105,  // quitte: 44 <-> mc: 44, class representative: 25
  105,  // quitte: 45 <-> mc: 45, class representative: 30
  46,   // quitte: 46 <-> mc: 46, class representative: 27
  47,   // quitte: 47 <-> mc: 47, class representative: 7
  48,   // quitte: 48 <-> mc: 48, class representative: 3
  49,   // quitte: 49 <-> mc: 49, class representative: 7
  50,   // quitte: 50 <-> mc: 50, class representative: 7
  51,   // quitte: 51 <-> mc: 51, class representative: 15
  105,  // quitte: 52 <-> mc: 56, class representative: 25
  53,   // quitte: 53 <-> mc: 57, class representative: 27
  105,  // quitte: 54 <-> mc: 58, class representative: 30
  55,   // quitte: 55 <-> mc: 59, class representative: 7
  105,  // quitte: 56 <-> mc: 52, class representative: 25
  105,  // quitte: 57 <-> mc: 53, class representative: 30
  58,   // quitte: 58 <-> mc: 54, class representative: 29
  59,   // quitte: 59 <-> mc: 55, class representative: 7
  105,  // quitte: 60 <-> mc: 60, class representative: 60
  105,  // quitte: 61 <-> mc: 61, class representative: 25
  105,  // quitte: 62 <-> mc: 62, class representative: 25
  63,   // quitte: 63 <-> mc: 63, class representative: 3
  64,   // quitte: 64 <-> mc: 128, class representative: 1
  105,  // quitte: 65 <-> mc: 129, class representative: 6
  105,  // quitte: 66 <-> mc: 130, class representative: 24
  105,  // quitte: 67 <-> mc: 131, class representative: 25
  68,   // quitte: 68 <-> mc: 136, class representative: 3
  69,   // quitte: 69 <-> mc: 137, class representative: 7
  105,  // quitte: 70 <-> mc: 138, class representative: 25
  71,   // quitte: 71 <-> mc: 139, class representative: 27
  105,  // quitte: 72 <-> mc: 132, class representative: 6
  105,  // quitte: 73 <-> mc: 133, class representative: 22
  105,  // quitte: 74 <-> mc: 134, class representative: 25
  105,  // quitte: 75 <-> mc: 135, class representative: 30
  76,   // quitte: 76 <-> mc: 140, class representative: 7
  77,   // quitte: 77 <-> mc: 141, class representative: 23
  78,   // quitte: 78 <-> mc: 142, class representative: 29
  79,   // quitte: 79 <-> mc: 143, class representative: 7
  80,   // quitte: 80 <-> mc: 144, class representative: 3
  81,   // quitte: 81 <-> mc: 145, class representative: 7
  105,  // quitte: 82 <-> mc: 146, class representative: 25
  83,   // quitte: 83 <-> mc: 147, class representative: 29
  84,   // quitte: 84 <-> mc: 152, class representative: 7
  85,   // quitte: 85 <-> mc: 153, class representative: 15
  105,  // quitte: 86 <-> mc: 154, class representative: 30
  87,   // quitte: 87 <-> mc: 155, class representative: 7
  105,  // quitte: 88 <-> mc: 148, class representative: 25
  105,  // quitte: 89 <-> mc: 149, class representative: 30
  105,  // quitte: 90 <-> mc: 150, class representative: 60
  105,  // quitte: 91 <-> mc: 151, class representative: 25
  92,   // quitte: 92 <-> mc: 156, class representative: 27
  93,   // quitte: 93 <-> mc: 157, class representative: 7
  105,  // quitte: 94 <-> mc: 158, class representative: 25
  95,   // quitte: 95 <-> mc: 159, class representative: 3
  105,  // quitte: 96 <-> mc: 160, class representative: 6
  105,  // quitte: 97 <-> mc: 161, class representative: 22
  105,  // quitte: 98 <-> mc: 162, class representative: 25
  105,  // quitte: 99 <-> mc: 163, class representative: 30
  105,  // quitte: 100 <-> mc: 168, class representative: 25
  105,  // quitte: 101 <-> mc: 169, class representative: 30
  105,  // quitte: 102 <-> mc: 170, class representative: 60
  105,  // quitte: 103 <-> mc: 171, class representative: 25
  105,  // quitte: 104 <-> mc: 164, class representative: 22
  105,  // quitte: 105 <-> mc: 165, class representative: 105
  105,  // quitte: 106 <-> mc: 166, class representative: 30
  105,  // quitte: 107 <-> mc: 167, class representative: 22
  105,  // quitte: 108 <-> mc: 172, class representative: 30
  105,  // quitte: 109 <-> mc: 173, class representative: 22
  105,  // quitte: 110 <-> mc: 174, class representative: 25
  105,  // quitte: 111 <-> mc: 175, class representative: 6
  112,  // quitte: 112 <-> mc: 176, class representative: 7
  113,  // quitte: 113 <-> mc: 177, class representative: 23
  114,  // quitte: 114 <-> mc: 178, class representative: 27
  115,  // quitte: 115 <-> mc: 179, class representative: 7
  116,  // quitte: 116 <-> mc: 184, class representative: 29
  117,  // quitte: 117 <-> mc: 185, class representative: 7
  105,  // quitte: 118 <-> mc: 186, class representative: 25
  119,  // quitte: 119 <-> mc: 187, class representative: 3
  105,  // quitte: 120 <-> mc: 180, class representative: 30
  105,  // quitte: 121 <-> mc: 181, class representative: 22
  105,  // quitte: 122 <-> mc: 182, class representative: 25
  105,  // quitte: 123 <-> mc: 183, class representative: 6
  105,  // quitte: 124 <-> mc: 188, class representative: 25
  105,  // quitte: 125 <-> mc: 189, class representative: 6
  105,  // quitte: 126 <-> mc: 190, class representative: 24
  127,  // quitte: 127 <-> mc: 191, class representative: 1
  128,  // quitte: 128 <-> mc: 64, class representative: 1
  105,  // quitte: 129 <-> mc: 65, class representative: 24
  105,  // quitte: 130 <-> mc: 66, class representative: 6
  105,  // quitte: 131 <-> mc: 67, class representative: 25
  105,  // quitte: 132 <-> mc: 72, class representative: 6
  105,  // quitte: 133 <-> mc: 73, class representative: 25
  105,  // quitte: 134 <-> mc: 74, class representative: 22
  105,  // quitte: 135 <-> mc: 75, class representative: 30
  136,  // quitte: 136 <-> mc: 68, class representative: 3
  105,  // quitte: 137 <-> mc: 69, class representative: 25
  138,  // quitte: 138 <-> mc: 70, class representative: 7
  139,  // quitte: 139 <-> mc: 71, class representative: 29
  140,  // quitte: 140 <-> mc: 76, class representative: 7
  141,  // quitte: 141 <-> mc: 77, class representative: 27
  142,  // quitte: 142 <-> mc: 78, class representative: 23
  143,  // quitte: 143 <-> mc: 79, class representative: 7
  105,  // quitte: 144 <-> mc: 80, class representative: 6
  105,  // quitte: 145 <-> mc: 81, class representative: 25
  105,  // quitte: 146 <-> mc: 82, class representative: 22
  105,  // quitte: 147 <-> mc: 83, class representative: 30
  105,  // quitte: 148 <-> mc: 88, class representative: 22
  105,  // quitte: 149 <-> mc: 89, class representative: 30
  105,  // quitte: 150 <-> mc: 90, class representative: 105
  105,  // quitte: 151 <-> mc: 91, class representative: 22
  105,  // quitte: 152 <-> mc: 84, class representative: 25
  105,  // quitte: 153 <-> mc: 85, class representative: 60
  105,  // quitte: 154 <-> mc: 86, class representative: 30
  105,  // quitte: 155 <-> mc: 87, class representative: 25
  105,  // quitte: 156 <-> mc: 92, class representative: 30
  105,  // quitte: 157 <-> mc: 93, class representative: 25
  105,  // quitte: 158 <-> mc: 94, class representative: 22
  105,  // quitte: 159 <-> mc: 95, class representative: 6
  160,  // quitte: 160 <-> mc: 96, class representative: 3
  105,  // quitte: 161 <-> mc: 97, class representative: 25
  162,  // quitte: 162 <-> mc: 98, class representative: 7
  163,  // quitte: 163 <-> mc: 99, class representative: 27
  105,  // quitte: 164 <-> mc: 104, class representative: 25
  105,  // quitte: 165 <-> mc: 105, class representative: 60
  105,  // quitte: 166 <-> mc: 106, class representative: 30
  105,  // quitte: 167 <-> mc: 107, class representative: 25
  168,  // quitte: 168 <-> mc: 100, class representative: 7
  105,  // quitte: 169 <-> mc: 101, class representative: 30
  170,  // quitte: 170 <-> mc: 102, class representative: 15
  171,  // quitte: 171 <-> mc: 103, class representative: 7
  172,  // quitte: 172 <-> mc: 108, class representative: 29
  105,  // quitte: 173 <-> mc: 109, class representative: 25
  174,  // quitte: 174 <-> mc: 110, class representative: 7
  175,  // quitte: 175 <-> mc: 111, class representative: 3
  176,  // quitte: 176 <-> mc: 112, class representative: 7
  177,  // quitte: 177 <-> mc: 113, class representative: 29
  178,  // quitte: 178 <-> mc: 114, class representative: 23
  179,  // quitte: 179 <-> mc: 115, class representative: 7
  105,  // quitte: 180 <-> mc: 120, class representative: 30
  105,  // quitte: 181 <-> mc: 121, class representative: 25
  105,  // quitte: 182 <-> mc: 122, class representative: 22
  105,  // quitte: 183 <-> mc: 123, class representative: 6
  184,  // quitte: 184 <-> mc: 116, class representative: 27
  105,  // quitte: 185 <-> mc: 117, class representative: 25
  186,  // quitte: 186 <-> mc: 118, class representative: 7
  187,  // quitte: 187 <-> mc: 119, class representative: 3
  105,  // quitte: 188 <-> mc: 124, class representative: 25
  105,  // quitte: 189 <-> mc: 125, class representative: 24
  105,  // quitte: 190 <-> mc: 126, class representative: 6
  191,  // quitte: 191 <-> mc: 127, class representative: 1
  192,  // quitte: 192 <-> mc: 192, class representative: 3
  105,  // quitte: 193 <-> mc: 193, class representative: 25
  105,  // quitte: 194 <-> mc: 194, class representative: 25
  105,  // quitte: 195 <-> mc: 195, class representative: 60
  196,  // quitte: 196 <-> mc: 200, class representative: 7
  197,  // quitte: 197 <-> mc: 201, class representative: 29
  105,  // quitte: 198 <-> mc: 202, class representative: 30
  105,  // quitte: 199 <-> mc: 203, class representative: 25
  200,  // quitte: 200 <-> mc: 196, class representative: 7
  105,  // quitte: 201 <-> mc: 197, class representative: 30
  202,  // quitte: 202 <-> mc: 198, class representative: 27
  105,  // quitte: 203 <-> mc: 199, class representative: 25
  204,  // quitte: 204 <-> mc: 204, class representative: 15
  205,  // quitte: 205 <-> mc: 205, class representative: 7
  206,  // quitte: 206 <-> mc: 206, class representative: 7
  207,  // quitte: 207 <-> mc: 207, class representative: 3
  208,  // quitte: 208 <-> mc: 208, class representative: 7
  209,  // quitte: 209 <-> mc: 209, class representative: 27
  105,  // quitte: 210 <-> mc: 210, class representative: 30
  105,  // quitte: 211 <-> mc: 211, class representative: 25
  212,  // quitte: 212 <-> mc: 216, class representative: 23
  213,  // quitte: 213 <-> mc: 217, class representative: 7
  105,  // quitte: 214 <-> mc: 218, class representative: 22
  105,  // quitte: 215 <-> mc: 219, class representative: 6
  216,  // quitte: 216 <-> mc: 212, class representative: 29
  105,  // quitte: 217 <-> mc: 213, class representative: 25
  105,  // quitte: 218 <-> mc: 214, class representative: 25
  105,  // quitte: 219 <-> mc: 215, class representative: 24
  220,  // quitte: 220 <-> mc: 220, class representative: 7
  221,  // quitte: 221 <-> mc: 221, class representative: 3
  105,  // quitte: 222 <-> mc: 222, class representative: 6
  223,  // quitte: 223 <-> mc: 223, class representative: 1
  224,  // quitte: 224 <-> mc: 224, class representative: 7
  105,  // quitte: 225 <-> mc: 225, class representative: 30
  226,  // quitte: 226 <-> mc: 226, class representative: 29
  105,  // quitte: 227 <-> mc: 227, class representative: 25
  228,  // quitte: 228 <-> mc: 232, class representative: 27
  105,  // quitte: 229 <-> mc: 233, class representative: 25
  105,  // quitte: 230 <-> mc: 234, class representative: 25
  105,  // quitte: 231 <-> mc: 235, class representative: 24
  232,  // quitte: 232 <-> mc: 228, class representative: 23
  105,  // quitte: 233 <-> mc: 229, class representative: 22
  234,  // quitte: 234 <-> mc: 230, class representative: 7
  105,  // quitte: 235 <-> mc: 231, class representative: 6
  236,  // quitte: 236 <-> mc: 236, class representative: 7
  105,  // quitte: 237 <-> mc: 237, class representative: 6
  238,  // quitte: 238 <-> mc: 238, class representative: 3
  239,  // quitte: 239 <-> mc: 239, class representative: 1
  240,  // quitte: 240 <-> mc: 240, class representative: 15
  241,  // quitte: 241 <-> mc: 241, class representative: 7
  242,  // quitte: 242 <-> mc: 242, class representative: 7
  243,  // quitte: 243 <-> mc: 243, class representative: 3
  244,  // quitte: 244 <-> mc: 248, class representative: 7
  245,  // quitte: 245 <-> mc: 249, class representative: 3
  105,  // quitte: 246 <-> mc: 250, class representative: 6
  247,  // quitte: 247 <-> mc: 251, class representative: 1
  248,  // quitte: 248 <-> mc: 244, class representative: 7
  105,  // quitte: 249 <-> mc: 245, class representative: 6
  250,  // quitte: 250 <-> mc: 246, class representative: 3
  251,  // quitte: 251 <-> mc: 247, class representative: 1
  252,  // quitte: 252 <-> mc: 252, class representative: 3
  253,  // quitte: 253 <-> mc: 253, class representative: 1
  254,  // quitte: 254 <-> mc: 254, class representative: 1
  255   // quitte: 255 <-> mc: 255, class representative: 0
};

} // namespace Cube_table
} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_TABLES_H