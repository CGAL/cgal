#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/Color.h>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <vector>

int main()
{
  using K = CGAL::Simple_cartesian<double>;
  using Color = CGAL::IO::Color;
  const std::filesystem::path filename("test_pmp_read_colors.off");
  std::ofstream output(filename);
  output << "COFF\n5 2 0\n"
            "0 0 0 255 0 0 255\n"
            "1 0 0 0 255 0 255\n"
            "1 1 0 0 0 255 255\n"
            "0 1 0 255 255 0 255\n"
            "0 0 0 255 0 255 255\n"
            "3 0 1 2 255 128 0 255\n"
            "3 4 2 3 0 128 255 255\n";
  output.close();

  const std::vector<Color> expected_vertices = {
    Color(255, 0, 0), Color(0, 255, 0), Color(0, 0, 255),
    Color(255, 255, 0), Color(255, 0, 255)
  };
  const std::vector<Color> expected_faces = {Color(255, 128, 0), Color(0, 128, 255)};
  for(bool repair : {false, true})
  {
    CGAL::Surface_mesh<K::Point_3> mesh;
    std::vector<Color> vertex_colors, face_colors;
    const bool read = CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, mesh,
      CGAL::parameters::repair_polygon_soup(repair)
        .vertex_color_output_iterator(std::back_inserter(vertex_colors))
        .face_color_output_iterator(std::back_inserter(face_colors)));
    assert(read);
    assert(vertex_colors == expected_vertices);
    assert(face_colors == expected_faces);
    assert(num_faces(mesh) == 2);
    if(repair)
      assert(num_vertices(mesh) == 4);
  }
  std::filesystem::remove(filename);
}
