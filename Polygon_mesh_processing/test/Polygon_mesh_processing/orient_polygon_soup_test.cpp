#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/algorithm.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
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
      std::size_t id = 0;
      scanner.scan_facet_vertex_index(id, j+1, i);
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

    CGAL::cpp98::random_shuffle(indices.begin(), indices.end());

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


template <class K, class Tag>
int test_pipeline()
{
  typedef typename K::Point_3 Point_3;
  typedef CGAL::Polyhedron_3<K> Polyhedron;

  std::vector<Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;
  Polyhedron ref1;

  shuffle_off("data/elephant.off", "elephant-shuffled.off");
  std::ifstream input("elephant-shuffled.off");
  if ( !input || !read_soup<K>(input, points, polygons)){
    std::cerr << "Error: can not shuffled file.\n";
    return 1;
  }
  input.close();
  input.open("data/elephant.off");
  if ( !input || !(input >> ref1)){
    std::cerr << "Error: can not read reference file.\n";
    return 1;
  }
  input.close();
  CGAL::Polygon_mesh_processing::orient_triangle_soup_with_reference_triangle_mesh<Tag>(ref1, points, polygons);

  CGAL::Polygon_mesh_processing::duplicate_non_manifold_edges_in_polygon_soup(points, polygons);

  Polyhedron poly;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(
        points, polygons, poly);
  typedef typename boost::property_map<Polyhedron, CGAL::dynamic_face_property_t<std::size_t> >::type Fccmap;
  Fccmap fim = get(CGAL::dynamic_face_property_t<std::size_t>(), poly);
  std::size_t id =0;
  for(auto f : faces(poly))
  {
    put(fim, f, id++);
  }
  CGAL::Polygon_mesh_processing::
      merge_reversible_connected_components(poly,
                                            CGAL::parameters::face_index_map(fim));

  Fccmap fccmap = get(CGAL::dynamic_face_property_t<std::size_t>(), poly);

  assert(CGAL::Polygon_mesh_processing::
         connected_components(poly, fccmap,
                              CGAL::parameters::face_index_map(fim)) == 1);
  return 0;
}

int main()
{
  assert(test_orient<Epic>(false) == 0);
  assert(test_orient<Epec>(false) == 0);

  int res = test_pipeline<Epic, CGAL::Sequential_tag>();
  assert(res == 0);
  res = test_pipeline<Epec, CGAL::Sequential_tag>();
  assert(res == 0);
#if defined(CGAL_LINKED_WITH_TBB)
  res = test_pipeline<Epic, CGAL::Parallel_tag>();
  assert(res == 0);
  //res = test_pipeline<Epec, CGAL::Parallel_tag>();
  //assert(res == 0);
#endif
  return 0;
}
