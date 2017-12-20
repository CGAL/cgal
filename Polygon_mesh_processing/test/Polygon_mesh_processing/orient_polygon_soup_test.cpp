#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Timer.h>

#include <fstream>
#include <algorithm>
#include <cstdlib>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Exact_predicates_exact_constructions_kernel Epec;

template <typename K>
std::istream& read_soup(
  std::istream& stream,
    std::vector<typename K::Point_3>& points,
  std::vector< std::vector<std::size_t> >& polygons)
{
  typedef typename K::Point_3 Point_3;
  CGAL::File_scanner_OFF scanner(stream);
  points.resize(scanner.size_of_vertices());
  polygons.resize(scanner.size_of_facets());
  for (std::size_t i = 0; i < scanner.size_of_vertices(); ++i)
  {
    double x, y, z, w;
    scanner.scan_vertex( x, y, z, w);
    points[i] = Point_3(x, y, z, w);
    scanner.skip_to_next_vertex( i);
  }
  if(!stream) { return stream; }

  for (std::size_t i = 0; i < scanner.size_of_facets(); ++i)
  {
    std::size_t no;
    scanner.scan_facet( no, i);
    polygons[i].resize(no);

    for(std::size_t j = 0; j < no; ++j) {
      std::size_t id;
      scanner.scan_facet_vertex_index(id, i);
      if(id < scanner.size_of_vertices())
      {
        polygons[i][j] = id;
      }
      else { return stream; }
    }
  }
  return stream;
}

void shuffle_off(const char* fname_in, const char* fname_out)
{
  std::ifstream input(fname_in);
  if ( !input ){
    std::cerr << "Error: can not read input file.\n";
    exit(1);
  }

  std::string OFF;
  int v, f, e;
  input >> OFF >> v >> f >> e;

  std::ofstream output(fname_out);
  output << "OFF\n" << v << " " << f << " 0\n\n";
  for(int i=0; i<v; ++i)
  {
    std::string line;
    while(line.empty())
      std::getline(input, line);
    output << line << "\n";
  }

  for (int i=0; i<f;++i)
  {
    int n;
    input >> n;
    std::vector<int> indices(n);
    for (int k=0;k<n;++k)
      input >> indices[k];

    std::random_shuffle(indices.begin(), indices.end());

    output << n;
    for (int k=0;k<n;++k)
      output << " " << indices[k];
    output << "\n";
  }
}

template <typename K>
int test_orient(const bool save_oriented) {
  typedef typename K::Point_3 Point_3;
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

  std::vector<Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;

  shuffle_off("data/elephant.off", "elephant-shuffled.off");
  std::ifstream input("elephant-shuffled.off");
  if ( !input || !read_soup<K>(input, points, polygons)){
    std::cerr << "Error: can not shuffled file.\n";
    return 1;
  }

  bool oriented
    = CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
  std::cerr << (oriented ? "Oriented." : "Not orientabled.") << std::endl;
  assert(oriented);

  if(oriented) {
    Surface_mesh mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(
      points, polygons, mesh);

    Polyhedron poly;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(
      points, polygons, poly);

    if (save_oriented)
    {
      std::ofstream out("elephant-oriented.off");
      out << poly;
      out.close();
    }
  }
  return 0;
}

int main()
{
  assert(test_orient<Epic>(false) == 0);
  assert(test_orient<Epec>(false) == 0);
  return 0;
}
