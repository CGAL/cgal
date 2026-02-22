// Copyright (c) 2025  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Rajdeep Singh

#include <CGAL/IO/MSH.h>
#include <CGAL/Simple_cartesian.h>

#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdio>

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------
typedef CGAL::Simple_cartesian<double>         Kernel;
typedef Kernel::Point_3                        Point;
typedef std::vector<Point>                     PointRange;
typedef std::vector<std::size_t>               Polygon;
typedef std::vector<Polygon>                   PolygonRange;

// ---------------------------------------------------------------------------
// Helper: build a minimal MSH 2.2 string in memory
// ---------------------------------------------------------------------------
static const std::string TRIANGLE_MESH_MSH = R"($MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 0.5 1.0 0.0
4 1.5 1.0 0.0
$EndNodes
$Elements
2
1 2 2 0 0 1 2 3
2 2 2 0 0 2 4 3
$EndElements
)";

// A mesh that contains a quad
static const std::string QUAD_MESH_MSH = R"($MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
4 0.0 1.0 0.0
$EndNodes
$Elements
1
1 3 2 0 0 1 2 3 4
$EndElements
)";

// A mesh with both surface triangles and volume tetrahedra (tets should be skipped)
static const std::string MIXED_MESH_MSH = R"($MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 0.5 1.0 0.0
4 0.5 0.5 1.0
$EndNodes
$Elements
2
1 2 2 0 0 1 2 3
2 4 2 0 0 1 2 3 4
$EndElements
)";

// A valid MSH with only non-surface elements (lines + points) — expect empty polygons
static const std::string LINES_ONLY_MSH = R"($MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
2
1 0.0 0.0 0.0
2 1.0 0.0 0.0
$EndNodes
$Elements
1
1 1 2 0 0 1 2
$EndElements
)";

// A completely empty file
static const std::string EMPTY_MSH = "";

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

void test_triangle_round_trip()
{
  std::cout << "test_triangle_round_trip ... ";

  PointRange   pts_in, pts_out;
  PolygonRange polys_in, polys_out;

  // Build a simple tetrahedron surface (4 triangles)
  pts_in = {
    Point(0, 0, 0),
    Point(1, 0, 0),
    Point(0, 1, 0),
    Point(0, 0, 1)
  };
  polys_in = {
    {0, 1, 2},
    {0, 1, 3},
    {0, 2, 3},
    {1, 2, 3}
  };

  // write to stream
  std::ostringstream oss;
  bool ok = CGAL::IO::write_MSH(oss, pts_in, polys_in);
  assert(ok && "write_MSH failed");

  // read back
  std::istringstream iss(oss.str());
  ok = CGAL::IO::read_MSH(iss, pts_out, polys_out);
  assert(ok && "read_MSH failed on round-trip data");
  assert(pts_out.size()   == pts_in.size()   && "point count mismatch");
  assert(polys_out.size() == polys_in.size() && "polygon count mismatch");

  // verify all faces have 3 nodes
  for(const auto& f : polys_out)
    assert(f.size() == 3);

  std::cout << "PASS\n";
}

void test_quad_round_trip()
{
  std::cout << "test_quad_round_trip ... ";

  PointRange   pts_in, pts_out;
  PolygonRange polys_in, polys_out;

  pts_in = {
    Point(0, 0, 0),
    Point(1, 0, 0),
    Point(1, 1, 0),
    Point(0, 1, 0)
  };
  polys_in = {
    {0, 1, 2, 3}   // one quad face
  };

  std::ostringstream oss;
  bool ok = CGAL::IO::write_MSH(oss, pts_in, polys_in);
  assert(ok && "write_MSH failed");

  std::istringstream iss(oss.str());
  ok = CGAL::IO::read_MSH(iss, pts_out, polys_out);
  assert(ok && "read_MSH failed");
  assert(pts_out.size()   == 4 && "point count");
  assert(polys_out.size() == 1 && "polygon count");
  assert(polys_out[0].size() == 4 && "quad size");

  std::cout << "PASS\n";
}

void test_read_triangle_mesh()
{
  std::cout << "test_read_triangle_mesh ... ";

  PointRange   pts;
  PolygonRange polys;

  std::istringstream iss(TRIANGLE_MESH_MSH);
  bool ok = CGAL::IO::read_MSH(iss, pts, polys);
  assert(ok);
  assert(pts.size()   == 4);
  assert(polys.size() == 2);
  for(const auto& f : polys)
    assert(f.size() == 3);

  std::cout << "PASS\n";
}

void test_read_quad_mesh()
{
  std::cout << "test_read_quad_mesh ... ";

  PointRange   pts;
  PolygonRange polys;

  std::istringstream iss(QUAD_MESH_MSH);
  bool ok = CGAL::IO::read_MSH(iss, pts, polys);
  assert(ok);
  assert(pts.size()   == 4);
  assert(polys.size() == 1);
  assert(polys[0].size() == 4);

  std::cout << "PASS\n";
}

void test_mixed_volume_surface()
{
  std::cout << "test_mixed_volume_surface ... ";

  PointRange   pts;
  PolygonRange polys;

  std::istringstream iss(MIXED_MESH_MSH);
  bool ok = CGAL::IO::read_MSH(iss, pts, polys);
  assert(ok);
  assert(pts.size()   == 4);
  // Only the triangle (type 2) should have been extracted; the tet (type 4) skipped
  assert(polys.size() == 1 && "volume element not skipped");
  assert(polys[0].size() == 3);

  std::cout << "PASS\n";
}

void test_lines_only()
{
  std::cout << "test_lines_only ... ";

  PointRange   pts;
  PolygonRange polys;

  std::istringstream iss(LINES_ONLY_MSH);
  // read_MSH should succeed (file is valid) but polygons stays empty
  bool ok = CGAL::IO::read_MSH(iss, pts, polys);
  assert(ok);
  assert(pts.size()   == 2);
  assert(polys.empty());

  std::cout << "PASS\n";
}

void test_empty_stream()
{
  std::cout << "test_empty_stream ... ";

  PointRange   pts;
  PolygonRange polys;

  std::istringstream iss(EMPTY_MSH);
  bool ok = CGAL::IO::read_MSH(iss, pts, polys);
  // An empty file has no $MeshFormat → should return false
  assert(!ok);

  std::cout << "PASS\n";
}

void test_nonexistent_file()
{
  std::cout << "test_nonexistent_file ... ";

  PointRange   pts;
  PolygonRange polys;

  bool ok = CGAL::IO::read_MSH(std::string("this_file_does_not_exist_12345678.msh"),
                                pts, polys,
                                CGAL::parameters::verbose(false));
  assert(!ok);

  std::cout << "PASS\n";
}

void test_file_overloads()
{
  std::cout << "test_file_overloads ... ";

  const std::string tmp = "cgal_test_msh_temp.msh";

  PointRange   pts_in = {
    Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1)
  };
  PolygonRange polys_in = {
    {0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}
  };

  bool ok = CGAL::IO::write_MSH(tmp, pts_in, polys_in);
  assert(ok && "file write failed");

  PointRange   pts_out;
  PolygonRange polys_out;
  ok = CGAL::IO::read_MSH(tmp, pts_out, polys_out);
  assert(ok && "file read failed");
  assert(pts_out.size()   == pts_in.size());
  assert(polys_out.size() == polys_in.size());

  std::remove(tmp.c_str());
  std::cout << "PASS\n";
}

void test_append_semantics()
{
  std::cout << "test_append_semantics ... ";

  PointRange   pts;
  PolygonRange polys;

  // First read
  {
    std::istringstream iss(TRIANGLE_MESH_MSH);
    bool ok = CGAL::IO::read_MSH(iss, pts, polys);
    assert(ok);
    assert(pts.size()   == 4);
    assert(polys.size() == 2);
  }

  // Second read — should append (not clear)
  {
    std::istringstream iss(TRIANGLE_MESH_MSH);
    bool ok = CGAL::IO::read_MSH(iss, pts, polys);
    assert(ok);
    assert(pts.size()   == 8 && "append: points should double");
    assert(polys.size() == 4 && "append: polygons should double");
  }

  // The second read's face indices should reference the appended points (index >= 4)
  for(std::size_t i = 2; i < 4; ++i)
    for(std::size_t j = 0; j < polys[i].size(); ++j)
      assert(polys[i][j] >= 4 && "appended face indices must reference the second set of nodes");

  std::cout << "PASS\n";
}

void test_fan_triangulation_on_write()
{
  std::cout << "test_fan_triangulation_on_write ... ";

  PointRange pts = {
    Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0),
    Point(0.5, 1.5, 0), Point(0, 1, 0)
  };
  PolygonRange polys = { {0, 1, 2, 3, 4} };  // pentagon → 3 triangles after fan-triangulation

  std::ostringstream oss;
  bool ok = CGAL::IO::write_MSH(oss, pts, polys);
  assert(ok);

  PointRange   pts_out;
  PolygonRange polys_out;
  std::istringstream iss(oss.str());
  ok = CGAL::IO::read_MSH(iss, pts_out, polys_out);
  assert(ok);
  assert(pts_out.size()   == 5 && "pentagon vertices");
  assert(polys_out.size() == 3 && "pentagon → 3 triangles");

  std::cout << "PASS\n";
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main()
{
  test_read_triangle_mesh();
  test_read_quad_mesh();
  test_mixed_volume_surface();
  test_lines_only();
  test_empty_stream();
  test_nonexistent_file();
  test_triangle_round_trip();
  test_quad_round_trip();
  test_file_overloads();
  test_append_semantics();
  test_fan_triangulation_on_write();

  std::cout << "\nAll MSH tests PASSED.\n";
  return 0;
}
