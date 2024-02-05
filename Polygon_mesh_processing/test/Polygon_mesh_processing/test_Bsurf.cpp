#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Bsurf/locally_shortest_path.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>


namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;
typedef PMP::Face_location<Mesh, double>                      Face_location;
typedef PMP::Edge_location<Mesh, double>                      Edge_location;


int main(int argc, char** argv)
{
  std::string filename = (argc > 1) ? std::string(argv[1])
    : CGAL::data_file_path("meshes/elephant.off");

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }
  PMP::Dual_geodesic_solver<double> solver;
  CGAL::Polygon_mesh_processing::init_geodesic_dual_solver(solver, mesh);


  std::size_t nb_hedges = halfedges(mesh).size();

  // take two random faces and pick the centroid
  CGAL::Random rnd = CGAL::get_default_random();
  //CGAL::Random rnd(1706604646);
  // CGAL::Random rnd(1695724381);
  // CGAL::Random rnd(1695813638);

  std::cout << "seed " << rnd.get_seed() << std::endl;
  Mesh::Halfedge_index h = *std::next(halfedges(mesh).begin(), rnd.get_int(0, nb_hedges));
  // Mesh::Halfedge_index h = *std::next(halfedges(mesh).begin(), 2*6178); // <---- interesting two different locally shortest paths

  std::cout << "h = " << h << "\n";

  std::ofstream out("locally_shortest_path.polylines.txt");

// test src/tgt being 2 opposite vertices of an edge
#if 1
  for (int loop=0; loop<3; ++loop)
  {
    Mesh::Vertex_index v1 = target(next(h, mesh), mesh);
    Mesh::Vertex_index v2 = target(next(opposite(h, mesh), mesh), mesh);

    bool first_run = true;
    std::size_t expected_size=0;
    for(Mesh::Halfedge_index h1 : CGAL::halfedges_around_target(v1, mesh))
    {
      Mesh::Face_index f1 = face(h1, mesh);
      Mesh::Halfedge_index hf1 = prev(halfedge(f1, mesh), mesh);
      int i1=0;
      while (target(hf1, mesh)!= v1)
      {
        hf1=next(hf1,mesh);
        ++i1;
      }
      Face_location src(f1, CGAL::make_array(0.,0.,0.));
      src.second[i1]=1;
      for(Mesh::Halfedge_index h2 : CGAL::halfedges_around_target(v2, mesh))
      {
        Mesh::Face_index f2 = face(h2, mesh);
        Mesh::Halfedge_index hf2 = prev(halfedge(f2, mesh), mesh);
        int i2=0;
        while (target(hf2, mesh)!= v2)
        {
          hf2=next(hf2,mesh);
          ++i2;
        }
        Face_location tgt(f2, CGAL::make_array(0.,0.,0.));
        tgt.second[i2]=1;

        std::cout << "Running " << f1 << " " << v1 << " | " << f2 << " " << v2 << "\n";

        std::vector<Edge_location> edge_locations;
        auto src_bk=src, tgt_bk=tgt;
        PMP::locally_shortest_path<double>(src, tgt, mesh, edge_locations, solver);
        assert(get<Mesh::Vertex_index>(PMP::get_descriptor_from_location(src,mesh))==v1);
        assert(get<Mesh::Vertex_index>(PMP::get_descriptor_from_location(tgt,mesh))==v2);

        out << edge_locations.size()+2;
        out << " " << PMP::construct_point(src, mesh);
        for (auto el : edge_locations)
          out << " " << PMP::construct_point(el, mesh);
        out << " " << PMP::construct_point(tgt, mesh) << "\n";
        out << std::flush;

        if (first_run)
        {
          first_run=false;
          expected_size=edge_locations.size();
        }
        if(edge_locations.size() != expected_size)
        {
          std::cout << edge_locations.size() << " vs " <<  expected_size << "\n";
        }
        CGAL_warning(edge_locations.size() == expected_size);

        src=src_bk;
        tgt=tgt_bk;
        edge_locations.clear();
        PMP::locally_shortest_path<double>(tgt, src, mesh, edge_locations, solver);
        assert(get<Mesh::Vertex_index>(PMP::get_descriptor_from_location(src,mesh))==v1);
        assert(get<Mesh::Vertex_index>(PMP::get_descriptor_from_location(tgt,mesh))==v2);

        out << edge_locations.size()+2;
        out << " " << PMP::construct_point(tgt, mesh);
        for (auto el : edge_locations)
          out << " " << PMP::construct_point(el, mesh);
        out << " " << PMP::construct_point(src, mesh) << "\n";
        out << std::flush;
        CGAL_warning(edge_locations.size() == expected_size);
      }
    }
    h = next(h, mesh);
  }
#endif

#if 1
  // test src is a vertex and tgt is on an edge
  for (int i=0; i<3; ++i)
  {
    Mesh::Vertex_index v1 = target(next(h, mesh), mesh);
    Mesh::Halfedge_index h2 = opposite(next(opposite(h, mesh), mesh), mesh);
    bool first_run=true;
    std::size_t expected_size=0;

    for(Mesh::Halfedge_index h1 : CGAL::halfedges_around_target(v1, mesh))
    {
      Mesh::Face_index f1 = face(h1, mesh);
      Mesh::Halfedge_index hf1 = prev(halfedge(f1, mesh), mesh);
      int i1=0;
      while (target(hf1, mesh)!= v1)
      {
        hf1=next(hf1,mesh);
        ++i1;
      }
      Face_location src(f1, CGAL::make_array(0.,0.,0.));
      src.second[i1]=1;

      // define tgt
      Mesh::Face_index f2 = face(h2, mesh);
      Mesh::Halfedge_index hf2 = prev(halfedge(f2, mesh), mesh);
      int k=0;
      while(hf2!=h2)
      {
        hf2=next(hf2, mesh);
        ++k;
      }

      Face_location tgt(f2, CGAL::make_array(0.5,0.5,0.5));
      tgt.second[(k+1)%3]=0;
      tgt.second[(k+2)%3]=0.25;
      tgt.second[(k)%3]=0.75;

      assert(get<Mesh::Vertex_index>(PMP::get_descriptor_from_location(src,mesh))==v1);
      assert(edge(get<Mesh::Halfedge_index>(PMP::get_descriptor_from_location(tgt,mesh)), mesh)==edge(h2,mesh));

      std::cout << "Running " << f1 << " " << v1 << " | " << f2 << " " << edge(h2,mesh) << "\n";
      // std::cout <<  "   " << PMP::construct_point(src, mesh) << " | " << PMP::construct_point(tgt, mesh) << "\n";
      std::vector<Edge_location> edge_locations;
      auto src_bk=src, tgt_bk=tgt;
      PMP::locally_shortest_path<double>(src, tgt, mesh, edge_locations, solver);

      out << edge_locations.size()+2;
      out << " " << PMP::construct_point(src, mesh);
      for (auto el : edge_locations)
        out << " " << PMP::construct_point(el, mesh);
      out << " " << PMP::construct_point(tgt, mesh) << "\n";
      out << std::flush;

      if (first_run)
      {
        first_run=false;
        expected_size=edge_locations.size();
      }
      CGAL_warning(edge_locations.size()==expected_size);

      assert(get<Mesh::Vertex_index>(PMP::get_descriptor_from_location(src,mesh))==v1);
      assert(edge(get<Mesh::Halfedge_index>(PMP::get_descriptor_from_location(tgt,mesh)), mesh)==edge(h2,mesh));

      src=src_bk;
      tgt=tgt_bk;
      edge_locations.clear();
      PMP::locally_shortest_path<double>(tgt, src, mesh, edge_locations, solver);

      out << edge_locations.size()+2;
      out << " " << PMP::construct_point(tgt, mesh);
      for (auto el : edge_locations)
        out << " " << PMP::construct_point(el, mesh);
      out << " " << PMP::construct_point(src, mesh) << "\n";
      out << std::flush;
      CGAL_warning(edge_locations.size()==expected_size);

      assert(get<Mesh::Vertex_index>(PMP::get_descriptor_from_location(src,mesh))==v1);
      assert(edge(get<Mesh::Halfedge_index>(PMP::get_descriptor_from_location(tgt,mesh)), mesh)==edge(h2,mesh));
    }
    h = next(h, mesh);
  }
#endif
#if 1
  // test src is on an edge and tgt is on an edge
  for (int i=0; i<3; ++i)
  {
    Mesh::Halfedge_index h1 = opposite(next(h, mesh), mesh);
    Mesh::Halfedge_index h2 = opposite(next(opposite(h, mesh), mesh), mesh);

    // define src
    Mesh::Face_index f1 = face(h1, mesh);
    Mesh::Halfedge_index hf1 = prev(halfedge(f1, mesh), mesh);
    int k1=0;
    while(hf1!=h1)
    {
      hf1=next(hf1, mesh);
      ++k1;
    }

    Face_location src(f1, CGAL::make_array(0.5,0.5,0.5));
    src.second[(k1+1)%3]=0;
    src.second[(k1+2)%3]=0.25;
    src.second[(k1)%3]=0.75;
    // define tgt
    Mesh::Face_index f2 = face(h2, mesh);
    Mesh::Halfedge_index hf2 = prev(halfedge(f2, mesh), mesh);
    int k2=0;
    while(hf2!=h2)
    {
      hf2=next(hf2, mesh);
      ++k2;
    }

    Face_location tgt(f2, CGAL::make_array(0.5,0.5,0.5));
    tgt.second[(k2+1)%3]=0;
    tgt.second[(k2+2)%3]=0.25;
    tgt.second[(k2)%3]=0.75;

    assert(edge(get<Mesh::Halfedge_index>(PMP::get_descriptor_from_location(src,mesh)), mesh)==edge(h1,mesh));
    assert(edge(get<Mesh::Halfedge_index>(PMP::get_descriptor_from_location(tgt,mesh)), mesh)==edge(h2,mesh));

    std::cout << "Running " << f1 << " " << edge(h1,mesh) << " | " << f2 << " " << edge(h2,mesh) << "\n";

    std::vector<Edge_location> edge_locations;
    auto src_bk=src, tgt_bk=tgt;
    PMP::locally_shortest_path<double>(src, tgt, mesh, edge_locations, solver);

    out << edge_locations.size()+2;
    out << " " << PMP::construct_point(src, mesh);
    for (auto el : edge_locations)
      out << " " << PMP::construct_point(el, mesh);
    out << " " << PMP::construct_point(tgt, mesh) << "\n";
    out << std::flush;

    std::size_t expected_size = edge_locations.size();

    assert(edge(get<Mesh::Halfedge_index>(PMP::get_descriptor_from_location(src,mesh)), mesh)==edge(h1,mesh));
    assert(edge(get<Mesh::Halfedge_index>(PMP::get_descriptor_from_location(tgt,mesh)), mesh)==edge(h2,mesh));

    src=src_bk;
    tgt=tgt_bk;
    edge_locations.clear();
    PMP::locally_shortest_path<double>(tgt, src, mesh, edge_locations, solver);

    out << edge_locations.size()+2;
    out << " " << PMP::construct_point(tgt, mesh);
    for (auto el : edge_locations)
      out << " " << PMP::construct_point(el, mesh);
    out << " " << PMP::construct_point(src, mesh) << "\n";
    out << std::flush;

    CGAL_warning(edge_locations.size()==expected_size);

    assert(edge(get<Mesh::Halfedge_index>(PMP::get_descriptor_from_location(src,mesh)), mesh)==edge(h1,mesh));
    assert(edge(get<Mesh::Halfedge_index>(PMP::get_descriptor_from_location(tgt,mesh)), mesh)==edge(h2,mesh));

    h = next(h, mesh);
  }
  #endif

  return 0;
}
