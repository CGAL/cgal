#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/GLTF.h>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3                Point;
typedef std::vector<std::size_t>       Face;

// ── helpers ──────────────────────────────────────────────────────────────────

static bool points_equal(const Point& a, const Point& b)
{
  return a.x() == b.x() && a.y() == b.y() && a.z() == b.z();
}

// ── tests ─────────────────────────────────────────────────────────────────────

// Write a hand-crafted tetrahedron, read it back, check exact values.
static void test_write_and_read_back()
{
  std::vector<Point> pts_in = {
    Point(0, 0, 0),
    Point(1, 0, 0),
    Point(0, 1, 0),
    Point(0, 0, 1)
  };
  std::vector<Face> faces_in = {
    {0, 1, 2},
    {0, 1, 3},
    {0, 2, 3},
    {1, 2, 3}
  };

  bool ok = CGAL::IO::write_GLTF("tmp.gltf", pts_in, faces_in);
  assert(ok && "write_GLTF failed");

  std::vector<Point> pts_out;
  std::vector<Face>  faces_out;
  ok = CGAL::IO::read_GLTF("tmp.gltf", pts_out, faces_out);
  assert(ok && "read_GLTF of written file failed");

  assert(pts_out.size()   == pts_in.size()   && "vertex count mismatch");
  assert(faces_out.size() == faces_in.size() && "face count mismatch");

  for (std::size_t i = 0; i < pts_in.size(); ++i)
    assert(points_equal(pts_in[i], pts_out[i]) && "vertex coordinate mismatch");

  for (std::size_t i = 0; i < faces_in.size(); ++i)
    assert(faces_out[i] == faces_in[i] && "face index mismatch");

  std::cout << "test_write_and_read_back: passed\n";
}

// Read Cube.gltf, write it, read it back, check counts and coordinates match.
static void test_roundtrip_cube(const std::string& cube_file)
{
  std::vector<Point> pts_orig;
  std::vector<Face>  faces_orig;

  bool ok = CGAL::IO::read_GLTF(cube_file, pts_orig, faces_orig);
  assert(ok && "failed to read Cube.gltf");

  ok = CGAL::IO::write_GLTF("tmp_cube.gltf", pts_orig, faces_orig);
  assert(ok && "failed to write tmp_cube.gltf");

  std::vector<Point> pts_rt;
  std::vector<Face>  faces_rt;
  ok = CGAL::IO::read_GLTF("tmp_cube.gltf", pts_rt, faces_rt);
  assert(ok && "failed to read tmp_cube.gltf");

  assert(pts_rt.size()   == pts_orig.size()   && "vertex count mismatch after round-trip");
  assert(faces_rt.size() == faces_orig.size() && "face count mismatch after round-trip");

  for (std::size_t i = 0; i < pts_orig.size(); ++i)
    assert(points_equal(pts_orig[i], pts_rt[i]) && "vertex mismatch after round-trip");

  for (std::size_t i = 0; i < faces_orig.size(); ++i)
    assert(faces_rt[i] == faces_orig[i] && "face mismatch after round-trip");

  std::cout << "test_roundtrip_cube: passed"
            << "  (" << pts_orig.size() << " vertices, "
            << faces_orig.size() << " faces)\n";
}

// Passing an empty mesh should succeed and produce a readable (empty) file.
static void test_write_empty()
{
  std::vector<Point> pts;
  std::vector<Face>  faces;

  bool ok = CGAL::IO::write_GLTF("tmp_empty.gltf", pts, faces);
  assert(ok && "write_GLTF with empty input failed");

  std::vector<Point> pts_rt;
  std::vector<Face>  faces_rt;
  ok = CGAL::IO::read_GLTF("tmp_empty.gltf", pts_rt, faces_rt);
  assert(ok && "read_GLTF of empty file failed");
  assert(pts_rt.empty()   && "expected no vertices");
  assert(faces_rt.empty() && "expected no faces");

  std::cout << "test_write_empty: passed\n";
}

// Read Lantern.gltf (external .bin buffer, 3 meshes, FLOAT positions,
// UNSIGNED_SHORT indices) and verify the expected geometry counts.
//
// Expected totals from the file:
//   LanternPole_Body:    926 verts, 2616 indices →  872 faces
//   LanternPole_Chain:   756 verts, 3744 indices → 1248 faces
//   LanternPole_Lantern: 2463 verts, 9822 indices → 3274 faces
//   Total:               4145 verts,               5394 faces
static void test_read_lantern(const std::string& lantern_file)
{
  std::vector<Point> points;
  std::vector<Face>  faces;

  bool ok = CGAL::IO::read_GLTF(lantern_file, points, faces);
  assert(ok && "failed to read Lantern.gltf");

  assert(points.size() == 4145 && "unexpected vertex count for Lantern");
  assert(faces.size()  == 5394 && "unexpected face count for Lantern");

  // Every face index must be a valid vertex index.
  for (const auto& f : faces)
    for (std::size_t idx : f)
      assert(idx < points.size() && "face index out of range");

  // All vertices must lie within the known bounding box of the model.
  // The largest extents come from the Body mesh accessor (with a small
  // margin for float-to-double conversion).
  const double margin = 0.01;
  const double x_max = 7.74559927  + margin;
  const double y_max = 12.8321095  + margin;
  const double z_max = 2.31570983  + margin;
  for (const auto& p : points)
  {
    assert(p.x() >= -x_max && p.x() <= x_max && "vertex x out of bounds");
    assert(p.y() >= -y_max && p.y() <= y_max && "vertex y out of bounds");
    assert(p.z() >= -z_max && p.z() <= z_max && "vertex z out of bounds");
  }

  std::cout << "test_read_lantern: passed"
            << "  (" << points.size() << " vertices, "
            << faces.size() << " faces)\n";
}

// Read Lantern.gltf, write it, read it back, and verify counts and
// coordinates are preserved.  Positions are originally FLOAT in the file;
// the reader promotes them to double, so the round-trip through our double
// writer is exact for those values.
static void test_roundtrip_lantern(const std::string& lantern_file)
{
  std::vector<Point> pts_orig;
  std::vector<Face>  faces_orig;

  bool ok = CGAL::IO::read_GLTF(lantern_file, pts_orig, faces_orig);
  assert(ok && "failed to read Lantern.gltf");

  ok = CGAL::IO::write_GLTF("tmp_lantern.gltf", pts_orig, faces_orig);
  assert(ok && "failed to write tmp_lantern.gltf");

  std::vector<Point> pts_rt;
  std::vector<Face>  faces_rt;
  ok = CGAL::IO::read_GLTF("tmp_lantern.gltf", pts_rt, faces_rt);
  assert(ok && "failed to read tmp_lantern.gltf");

  assert(pts_rt.size()   == pts_orig.size()   && "vertex count mismatch after Lantern round-trip");
  assert(faces_rt.size() == faces_orig.size() && "face count mismatch after Lantern round-trip");

  for (std::size_t i = 0; i < pts_orig.size(); ++i)
    assert(points_equal(pts_orig[i], pts_rt[i]) && "vertex mismatch after Lantern round-trip");

  for (std::size_t i = 0; i < faces_orig.size(); ++i)
    assert(faces_rt[i] == faces_orig[i] && "face index mismatch after Lantern round-trip");

  std::cout << "test_roundtrip_lantern: passed"
            << "  (" << pts_orig.size() << " vertices, "
            << faces_orig.size() << " faces)\n";
}

// ── main ─────────────────────────────────────────────────────────────────────

int main(int argc, char** argv)
{
  const std::string cube_file    = (argc > 1) ? argv[1] : "data/Cube.gltf";
  const std::string lantern_file = (argc > 2) ? argv[2] : "data/Lantern.gltf";

  test_write_and_read_back();
  test_roundtrip_cube(cube_file);
  test_write_empty();
  test_read_lantern(lantern_file);
  test_roundtrip_lantern(lantern_file);

  std::cout << "All tests passed.\n";
  return EXIT_SUCCESS;
}
