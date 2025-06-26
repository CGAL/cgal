#include <CGAL/Simple_cartesian.h>
#include <CGAL/iterator.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/TDS_3/Indexed_storage.h>
#include <sstream>

int main() {
  using K = CGAL::Simple_cartesian<double>;
  using Weighted_point = K::Weighted_point_3;
  using Vb = CGAL::VertexWithPoint<K>;
  using Cb = CGAL::Cell4Regular<K>;
  using TDS = CGAL::Triangulation_data_structure_3<Vb, Cb,
                                                   CGAL::Sequential_tag,
                                                   CGAL::Index_tag>;
  using Vertex_handle = TDS::Vertex_handle;
  using Vertex_iterator = TDS::Vertex_iterator;
  using Cell_handle = TDS::Cell_handle;
  using Cell_circulator = TDS::Cell_circulator;
  using Facet_circulator = TDS::Facet_circulator;
  using Edge = TDS::Edge;

  TDS tds;

  Vertex_handle inf, v0, v1, v2, v3;
  inf = tds.insert_first_finite_cell(v0, v1, v2, v3);

  for(const auto& v : tds.vertices()) {
    std::cout << v << std::endl;
  }

  for(Vertex_iterator vit = tds.vertices_begin(); vit != tds.vertices_end(); ++vit) {
    std::cout << vit->point() << std::endl;
  }

  for(auto c : tds.cells()) {
    std::cout << c.index() << std::endl;
    c.hide_point(Weighted_point{});
  }

  std::vector<Cell_handle> incident_cells;
  tds.incident_cells_3(inf, incident_cells);
  assert(incident_cells.size() == 4);

  for(auto f : tds.facets()) {
    std::cout << f.first << " " << f.second << std::endl;
  }

  for(auto e : tds.edges()) {
    std::cout << e.first << " " << e.second << " " << e.third << std::endl;
  }

  Cell_handle c = inf->cell();
  Edge e(c, 0, 1);

  std::cout << "Incident cells to edge:\n";
  Cell_circulator cc = tds.incident_cells(e), cdone(cc);
  do {
    // cc->neighbor(0);
    Cell_handle ch = cc;
    std::cout << ch << " ";
    ++cc;
  } while (cc != cdone);
  std::cout << std::endl;

  std::cout << "Incident facets to edge:\n";
  Facet_circulator fc = tds.incident_facets(e), fdone(fc);
  do {
    Cell_handle ch = fc->first;
    std::cout << ch << " ";
    ++fc;
  } while (fc != fdone);
  std::cout << std::endl;


  struct IgnoreZero {
    bool operator()(const Vertex_handle& v) const {
      return v->index().id() == 0; // Example filter: skip vertex with index 0
    }
  };
  using Ignore_zero_iterator = CGAL::Filter_iterator<Vertex_iterator,IgnoreZero>;

  Ignore_zero_iterator vit(tds.vertices_begin(),IgnoreZero());
  ++vit;
  vit->index().id();

 // What is different from operator>> for Triangulation_3
  std::stringstream point("0 1 2");
  std::vector<Vertex_handle> vertices;
  vertices.push_back(v0);
  point >> *vertices[0];
  std::cout << v0->point() << std::endl;

  return 0;
}
