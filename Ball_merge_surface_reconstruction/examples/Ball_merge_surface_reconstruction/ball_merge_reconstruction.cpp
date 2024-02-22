#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Ball_merge_surface_reconstruction.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;



int main(int argc, char **argv)
{
  double tlen;
  const char *inFilename = argv[1];//Filename
  std::ifstream inStream(inFilename);//Read the file
  std::vector<std::array<unsigned char, 3>> meshVertexColors;//Variable to hold if the input has textures
  Point p;
  double par = atof(argv[2]);//Parameter to check IR
  int option = atoi(argv[3]);//Option to select global and local variants - 1 for global and 0 for local
  if (argc >= 5)//If local, there is an option to give extra optional parameter to filter long triangles
    tlen = atof(argv[4]);
  else
    tlen = 200.0;//If not provided, default 200 will be taken

  int vIndex = 0;
  std::vector<std::pair<Point, unsigned>> points;

  while (!inStream.eof()){//Read the data from the input file - coordinates, texture information (if available)
    std::string line;
    std::getline(inStream, line);
    std::istringstream str(line);
    if (str >> p){
      int r, g, b;
      if (str >> r){
        str >> g;
        str >> b;
        meshVertexColors.push_back({(unsigned char)r, (unsigned char)g, (unsigned char)b});
      }
      points.push_back(std::make_pair(Point(p.x(), p.y(), p.z()), vIndex++));
    }
  }

  CGAL::Ball_merge_surface_reconstruction<K, CGAL::Parallel_tag> bmsr;
  bmsr.option=option;

  bmsr(points, par, tlen);

  // AF: In case of duplicated points the colors of vertices will be wrong

  std::vector<Point> meshVertexPositions(points.size());//Preparing the tetrahedra for creating the PLY file & PLY file writing starts
  std::vector<std::vector<int>> meshFaceIndices;

  bmsr.result(meshVertexPositions, meshFaceIndices);

  std::string st = std::to_string(par) + "out"+std::to_string(tlen)+".ply";

  CGAL::IO::write_polygon_soup(st, meshVertexPositions, meshFaceIndices);
std::cout << "#faces " << meshFaceIndices.size() << "\n";
  meshFaceIndices.clear();

  if (option == 1){//Sometimes, in the gloabl case, the largest group would be a mould created by the convex hull, just to avoid it, we will write the second largest group as well
    st = st + "1.ply";

    for (auto vit = bmsr.dt3.finite_cells_begin(); vit != bmsr.dt3.finite_cells_end(); vit++)
      for (int i = 0; i < 4; i++)
        if (vit->info() == bmsr.secondgroup && (vit->neighbor(i)->info() != bmsr.secondgroup || bmsr.dt3.is_infinite(vit->neighbor(i)))){//Write the triangles between cells if they have have different labels and one of them is labeled as the same as the second largest group
          std::vector<int> indices(3);
          for (int j = 0; j < 3; j++)
            indices[j] = vit->vertex((i + 1 + j) % 4)->info();
          meshFaceIndices.push_back(indices);
        }
      CGAL::IO::write_polygon_soup(st, meshVertexPositions, meshFaceIndices);
  }

  return 0;
}
