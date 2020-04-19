#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#define CGAL_GENERIC_P2T2 // @todo still needed but to remove eventually

// #define CGAL_DEBUG_P2T2

#include <CGAL/internal/Generic_P2T2/Periodic_2_Delaunay_triangulation_2_generic.h>

#include <CGAL/point_generators_2.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                            K;
typedef K::FT                                                                          FT;
typedef K::Vector_2                                                                    Vector;

typedef typename CGAL::Periodic_2_offset_2                                             Offset;
typedef CGAL::Lattice_2<K>                                                             Lattice;
typedef CGAL::Periodic_2_triangulations_2::internal::Lattice_construct_point_2<K>      CP2;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_base_2<K, Offset, Lattice, CP2> GT;

typedef CGAL::Periodic_2_triangulation_vertex_base_2_generic<GT>                       Vb;
typedef CGAL::Triangulation_face_base_with_info_2<std::pair<bool, CGAL::Color>, GT>    Fbb;
typedef CGAL::Periodic_2_triangulation_face_base_2_generic<GT, Fbb>                    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                                   Tds;
typedef CGAL::Periodic_2_Delaunay_triangulation_2_generic<GT, Tds>                     GPDT;

typedef GPDT::Point                                                                    Point;

typedef GPDT::Vertex_handle                                                            Vertex_handle;
typedef GPDT::Face_handle                                                              Face_handle;

typedef CGAL::Creator_uniform_2<FT, Point>                        Creator;
typedef CGAL::Random_points_in_square_2<Point, Creator>           Point_generator;

// Remark: the basis
// basis = CGAL::make_array(Vector(4, 1), Vector(-2.5, -1));
// with the point set
// std::vector<Point> pts { Point(0, 0), Point(-0.4, -0.6) };
// are not in general position. If using inexact constructions,
// this yields an inconsistent triangulation (i.e. the are 4
// points forming a rectangle, but different periodic copies
// of this rectangle are triangulated differently), making
// the triangulation break down.
// --> Use exact constructions for now, in the future hopefully
// we replace these constructions with predicates.
void basic_gpdt()
{
  std::array<Vector, 2> basis;

  // the reduced version of this basis should be equivalent to
  // CGAL::make_array(Vector(-0.5, 1), Vector(1.5, 0));
  basis = CGAL::make_array(Vector(4, 1), Vector(-2.5, -1));

  std::vector<Point> pts { Point(0, 0), Point(-0.21, -0.6) };

  GPDT Tr(pts.begin(), pts.end(), basis);
  GPDT::Vertex_handle vh = Tr.insert(Point(12345.6, 5432.1));
  GPDT::Vertex_handle vh2 = Tr.insert(Point(0, 0.1));
  GPDT::Vertex_handle vh3 = Tr.insert(Point(0.5, -0.4));
  GPDT::Vertex_handle vh4 = Tr.insert(Point(0.2, -0.2));
  GPDT::Vertex_handle vh5 = Tr.insert(Point(0, -0.52));
  GPDT::Vertex_handle vh6 = Tr.insert(Point(0.6, 0));
  std::cout << "Is simplicial complex (vh6): " << Tr.is_simplicial_complex() << std::endl;
  std::cout << "Too big faces (vh6): " << Tr.too_big_faces() << std::endl;
  GPDT::Vertex_handle vh7 = Tr.insert(Point(0.8, 0.55));
  std::cout << "Is simplicial complex (vh7): " << Tr.is_simplicial_complex() << std::endl;
  std::cout << "Too big faces (vh7): " << Tr.too_big_faces() << std::endl;

  std::cout << "Number of vertices: " << Tr.number_of_vertices() << std::endl;
  std::cout << "Number of edges: " << Tr.number_of_edges() << std::endl;
  std::cout << "Number of faces: " << Tr.number_of_faces() << std::endl;

  for(auto ppit = Tr.periodic_points_begin(); ppit != Tr.periodic_points_end(); ++ppit)
    std::cout << ppit->first << "  " << ppit->second << std::endl;

  // Draw everything surrounding the vertex in red.
  GPDT::Face_handle some_face;
  std::set<GPDT::Face_handle> incident_faces2 = Tr.incident_faces(vh7);
  for(GPDT::Face_handle fh : incident_faces2)
  {
    fh->info().first = true;
    fh->info().second = CGAL::Color(255, 0, 0);
    some_face = fh;
  }
  some_face->info().first = true;
  some_face->info().second = CGAL::Color(255, 0, 192);

  // Make the neighbours of the cell blueish.
  for(int i=0; i<3; ++i)
  {
    GPDT::Face_handle nbfh = Tr.neighbor(some_face, i);

    if(nbfh->info().first)
    {
      const CGAL::Color& c = nbfh->info().second;
      nbfh->info().second = CGAL::Color(c.red(), c.green(), 128);
    }
    else
    {
      nbfh->info().first = true;
      nbfh->info().second = CGAL::Color(255, 255, 128);
    }
  }

  std::ofstream out_dt2("dt2.off");
  CGAL::draw_t2(out_dt2, Tr.dt2);
  out_dt2.close();

  Tr.convert_to_1_cover();

  for(auto vit=Tr.p2dt2.vertices_begin(); vit!=Tr.p2dt2.vertices_end(); ++vit)
    std::cout << vit->point() << " 0" << std::endl;

  std::ofstream out_b("before_p2dt2.off");
  CGAL::write_P2T2_to_OFF(out_b, Tr.p2dt2);
  out_b.close();

  std::cout << "Tr.p2t2.nv: " << Tr.p2dt2.number_of_vertices() << std::endl;
  std::cout << "Tr.p2t2.nf: " << Tr.p2dt2.number_of_faces() << std::endl;

  Tr.insert(Point(0.3, 0.12));

  std::cout << "Tr.p2t2.nv: " << Tr.p2dt2.number_of_vertices() << std::endl;
  std::cout << "Tr.p2t2.nf: " << Tr.p2dt2.number_of_faces() << std::endl;

  CGAL_assertion(Tr.p2dt2.is_valid());
  std::cout << "All good!" << std::endl;

  std::ofstream out_a("after_p2dt2.off");
  CGAL::write_P2T2_to_OFF(out_a, Tr.p2dt2);
}

std::vector<Point> generate_lattice()
{
  std::vector<Point> pts;

  for(int i=0; i<10; ++i)
    for(int j=0; j<10; ++j)
      pts.emplace_back(i * 0.14, j * 0.32);

  return pts;
}

std::vector<Point> generate_random_points(std::size_t n, CGAL::Random& rnd)
{
  std::vector<Point> pts;
  pts.reserve(n);

  std::copy_n(Point_generator(1, rnd), n, back_inserter(pts)); // square centered on 0, side 2

  return pts;
}

int main(int argc, char** argv)
{
  const std::size_t number_of_points = (argc > 1) ? std::atoi(argv[1]) : 800;
  std::cout << number_of_points << " random points" << std::endl;
  const int random_seed = (argc > 2) ? std::atoi(argv[2]) : CGAL::get_default_random().get_int(0, (1 << 30));
  std::cout << "random seed: " << random_seed << std::endl;

  CGAL::Random rnd(random_seed);
  const std::vector<Point> pts = generate_random_points(number_of_points, rnd);

  // the reduced version of this basis should be equivalent to
  // CGAL::make_array(Vector(-0.5, 1), Vector(1.5, 0));
  const std::array<Vector, 2> basis = CGAL::make_array(Vector(4, 1), Vector(-2.5, -1));
  const CGAL::Lattice_2<K> lattice(basis);

  GPDT Tr(lattice);
  for(const Point& pt : pts)
  {
    if(lattice.is_in_scaled_domain(pt))
      Tr.insert(pt);
  }

  std::cout << "Number of vertices: " << Tr.number_of_vertices() << std::endl;
  std::cout << "Number of edges: " << Tr.number_of_edges() << std::endl;
  std::cout << "Number of faces: " << Tr.number_of_faces() << std::endl;

  if(Tr.is_1_cover())
  {
    std::ofstream out("final.off");
    CGAL::write_P2T2_to_OFF(out, Tr.p2dt2);
  }
  else
  {
    std::ofstream out("final.off");
    CGAL::draw_t2(out, Tr.dt2);
  }

  return EXIT_SUCCESS;
}
