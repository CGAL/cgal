#include <CGAL/Simple_cartesian.h>
#include <CGAL/iterator.h>
#include <CGAL/TDS_3/Indexed_storage.h>

int main() {
  using K = CGAL::Simple_cartesian<double>;
  using Point = K::Point_3;
  using Vb = CGAL::VertexWithPoint<K>;
  using Cb = CGAL::Cell4Regular<K>;
  using TDS = CGAL::Indexed_storage<Vb, Cb>;
  TDS tds;
  using Vertex_handle = TDS::Vertex_handle;
  using Vertex_iterator = TDS::Vertex_iterator;
  using Cell_handle = TDS::Cell_handle;
  using Cell_circulator = TDS::Cell_circulator;
  using Facet_circulator = TDS::Facet_circulator;
  using Edge = TDS::Edge;

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
    c.hide_point(Point{});
  }

  inf->cell()->storage().my_other_vertex = v0;
  inf->cell()->storage().my_other_cell = v1->cell();

  std::vector<Cell_handle> incident_cells;
  tds.just_incident_cells_3(inf, incident_cells);
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


  return 0;
}
