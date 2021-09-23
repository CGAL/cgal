#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Random.h>

// Define the kernel.
typedef CGAL::Exact_predicates_exact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                       Point_3;
typedef Kernel::Vector_3                                      Vector_3;
typedef Kernel::Segment_3                                     Segment_3;

// Define the structure.
typedef CGAL::Delaunay_triangulation_3< Kernel >    DT;

typedef DT::Vertex_handle                           Vertex_handle;
typedef DT::Cell_handle                             Cell_handle;
typedef DT::Edge                                    Edge;
typedef DT::Facet                                   Facet;
typedef DT::Segment_simplex_iterator                Segment_simplex_iterator;

void test_vertex_edge_vertex(const DT& dt, const std::size_t& nb_tests)
{
  std::cout << "* test_vertex_edge_vertex *" << std::endl;
  std::vector<Edge> edges;
  for (DT::Finite_edges_iterator eit = dt.finite_edges_begin();
       eit != dt.finite_edges_end() && edges.size() < nb_tests;
       ++eit)
  {
    edges.push_back(*eit);
  }

  for (std::size_t i = 0; i < nb_tests; ++i)
  {
    Vertex_handle v1 = edges[i].first->vertex(edges[i].second);
    Vertex_handle v2 = edges[i].first->vertex(edges[i].third);
    Vector_3 v(v1->point(), v2->point());

    std::cout << "TEST " << i << " (" << v1->point()
                          << " ** " << v2->point() <<")"
                          << std::endl;
    std::cout << "\t(";
    Segment_simplex_iterator st
      = dt.segment_traverser_simplices_begin((v1->point() - 2.*v),
                                             (v2->point() + 3.*v));
    Segment_simplex_iterator stend
      = dt.segment_traverser_simplices_end();
    for (; st != stend; ++st)
    {
      std::cout << st->dimension();
      if(st->dimension() == 3
        && Cell_handle(*st) != Cell_handle()
        && dt.is_infinite(Cell_handle(*st)))
        std::cout << "i";
      std::cout << " ";

      if (st->dimension() == 0 && Vertex_handle(*st) == v1)
      {
        ++st;
        std::cout << st->dimension() << " ";
        assert(st->dimension() == 1);
        Edge e(*st);
        Vertex_handle ve1 = e.first->vertex(e.second);
        Vertex_handle ve2 = e.first->vertex(e.third);
        assert((ve1 == v1 && ve2 == v2)
            || (ve1 == v2 && ve2 == v1));

        ++st;
        std::cout << st->dimension() << " ";
        assert(st->dimension() == 0);
        assert(Vertex_handle(*st) == v2);
      }
    }
    std::cout << ")" << std::endl;
  }
}

void test_edge_facet_edge(const DT& dt, const std::size_t& nb_tests)
{
  std::cout << "* test_edge_facet_edge *" << std::endl;
  std::vector<Facet> facets;
  for (DT::Finite_facets_iterator fit = dt.finite_facets_begin();
    fit != dt.finite_facets_end() && facets.size() < nb_tests;
    ++fit)
  {
    facets.push_back(*fit);
  }
  for (std::size_t i = 0; i < nb_tests; ++i)
  {
    const int fi = facets[i].second;
    Vertex_handle v1 = facets[i].first->vertex((fi + 1) % 4);
    Vertex_handle v2 = facets[i].first->vertex((fi + 2) % 4);
    Vertex_handle v3 = facets[i].first->vertex((fi + 3) % 4);

    Point_3 p1 = CGAL::midpoint(v1->point(), v2->point());
    Point_3 p2 = CGAL::midpoint(v2->point(), v3->point());
    Vector_3 v(p1, p2);

    std::cout << "TEST " << i << " (" << p1 << " ** " << p2 << ")"
      << std::endl;
    std::cout << "\t(";
    Segment_simplex_iterator st
      = dt.segment_traverser_simplices_begin((p1 - 2. * v), (p2 + 3.  * v));
    Segment_simplex_iterator stend
      = dt.segment_traverser_simplices_end();
    for (; st != stend; ++st)
    {
      std::cout << st->dimension();
      if (st->dimension() == 3
        && Cell_handle(*st) != Cell_handle()
        && dt.is_infinite(Cell_handle(*st)))
        std::cout << "i";
      std::cout << " ";

      if (st->dimension() == 1)
      {
        Edge e = *st;
        Vertex_handle va = e.first->vertex(e.second);
        Vertex_handle vb = e.first->vertex(e.third);
        if ((va == v1 && vb == v2)
          || (va == v2 && vb == v1))
        {
          ++st;
          std::cout << st->dimension() << " ";
          assert(st->dimension() == 2);

          ++st;
          std::cout << st->dimension() << " ";
          assert(st->dimension() == 1);
          Edge e2 = *st;
          Vertex_handle va2 = e2.first->vertex(e2.second);
          Vertex_handle vb2 = e2.first->vertex(e2.third);
          assert(va == va2 || va == vb2 || vb == va2 || vb == vb2);
        }
      }
    }
    std::cout << ")" << std::endl;
  }
}

void test_edge_facet_vertex(const DT& dt, const std::size_t& nb_tests)
{
  std::cout << "* test_edge_facet_vertex *" << std::endl;
  std::vector<Facet> facets;
  for (DT::Finite_facets_iterator fit = dt.finite_facets_begin();
    fit != dt.finite_facets_end() && facets.size() < nb_tests;
    ++fit)
  {
    facets.push_back(*fit);
  }
  for (std::size_t i = 0; i < nb_tests; ++i)
  {
    const int fi = facets[i].second;
    Vertex_handle v1 = facets[i].first->vertex((fi + 1) % 4);
    Vertex_handle v2 = facets[i].first->vertex((fi + 2) % 4);
    Vertex_handle v3 = facets[i].first->vertex((fi + 3) % 4);

    Point_3 p1 = CGAL::midpoint(v1->point(), v2->point());
    Point_3 p2 = v3->point();
    Vector_3 v(p1, p2);

    std::cout << "TEST " << i << " (" << p1 << " ** " << p2 << ")"
      << std::endl;
    std::cout << "\t(";
    Segment_simplex_iterator st
      = dt.segment_traverser_simplices_begin((p1 - 2. * v), (p2 + 3. * v));
    Segment_simplex_iterator stend
      = dt.segment_traverser_simplices_end();
    for (; st != stend; ++st)
    {
      std::cout << st->dimension();
      if (st->dimension() == 3 && dt.is_infinite(Cell_handle(*st)))
        std::cout << "i";
      std::cout << " ";

      if (st->dimension() == 1)
      {
        Edge e = *st;
        Vertex_handle va = e.first->vertex(e.second);
        Vertex_handle vb = e.first->vertex(e.third);
        if ((va == v1 && vb == v2) || (va == v2 && vb == v1))
        {
          ++st;
          std::cout << st->dimension() << " ";
          assert(st->dimension() == 2);

          ++st;
          std::cout << st->dimension() << " ";
          assert(st->dimension() == 0);
          assert(Vertex_handle(*st) == v3);
        }
        ++st;
        std::cout << st->dimension() << " ";
        if (st == stend)
          break;
        else if (st->dimension() == 3
          && Cell_handle(*st) != Cell_handle()
          && dt.is_infinite(Cell_handle(*st)))
          std::cout << "i ";
        assert(st->dimension() == 3);
      }
    }
    std::cout << ")" << std::endl;
  }
}

void test_vertex_facet_edge(const DT& dt, const std::size_t& nb_tests)
{
  std::cout << "* test_vertex_facet_edge *" << std::endl;
  std::vector<Facet> facets;
  DT::Finite_facets_iterator fit = dt.finite_facets_begin();
  ++fit; ++fit; ++fit; //just avoid using the same faces as for test_edge_facet_vertex
  for (; fit != dt.finite_facets_end() && facets.size() < nb_tests;
       ++fit)
  {
    facets.push_back(*fit);
  }
  for (std::size_t i = 0; i < nb_tests; ++i)
  {
    const int fi = facets[i].second;
    Vertex_handle v1 = facets[i].first->vertex((fi + 1) % 4);
    Vertex_handle v2 = facets[i].first->vertex((fi + 2) % 4);
    Vertex_handle v3 = facets[i].first->vertex((fi + 3) % 4);

    Point_3 p1 = v1->point();
    Point_3 p2 = CGAL::midpoint(v2->point(), v3->point());
    Vector_3 v(p1, p2);

    std::cout << "TEST " << i << " (" << p1 << " ** " << p2 << ")"
      << std::endl;
    std::cout << "\t(";
    Segment_simplex_iterator st = dt.segment_traverser_simplices_begin(p1 - 2.*v, p2 + 3.*v);
    Segment_simplex_iterator stend = dt.segment_traverser_simplices_begin(p1 - 2. * v, p2 + 3. * v);
    for (; st != stend; ++st)
    {
      std::cout << st->dimension();
      if (st->dimension() == 3
        && Cell_handle(*st) != Cell_handle()
        && dt.is_infinite(Cell_handle(*st)))
        std::cout << "i";
      std::cout << " ";

      if (st->dimension() == 0 && Vertex_handle(*st) == v1)
      {
        ++st;
        std::cout << st->dimension() << " ";
        assert(st->dimension() == 2);
        assert(Facet(*st) == facets[i]
            || Facet(*st) == dt.mirror_facet(facets[i]));
        ++st;
        std::cout << st->dimension() << " ";
        assert(st->dimension() == 1);
        Edge e(*st);
        Vertex_handle va = e.first->vertex(e.second);
        Vertex_handle vb = e.first->vertex(e.third);
        assert((va == v2 && vb == v3) || (va == v3 && vb == v2));
      }
    }
    std::cout << ")" << std::endl;
  }
}

void test_triangulation_on_a_grid()
{
  std::cout << "* test_triangulation_on_a_grid *" << std::endl;
  DT dt;
  for (double x = 0.; x < 11.; x = x + 1.)
    for (double y = 0.; y < 11.; y = y + 1.)
      for (double z = 0.; z < 11.; z = z + 1.)
        dt.insert(Point_3(x, y, z));

  int nb_queries = 5;
  std::vector<Segment_3> queries(nb_queries);
  //along an axis of the grid
  queries[0] = Segment_3(Point_3(1., 1., 1.), Point_3(1., 8., 1.));
  //along an axis, but between two layers
  queries[1] = Segment_3(Point_3(1., 1.5, 1.), Point_3(10., 1.5, 1.));
  //along a diagonal
  queries[2] = Segment_3(Point_3(1., 1., 1.), Point_3(6., 6., 6.));
  //along a diagonal, in the plane (y = 1)
  queries[3] = Segment_3(Point_3(1., 1., 1.), Point_3(7., 1., 7.));
  //along a border of the cube
  queries[4] = Segment_3(Point_3(0., 0., 0.), Point_3(11., 0., 5.));

  for(const Segment_3& s : queries)
  {
    std::cout << "Query segment : (" << s.source()
                         << ") to (" << s.target() << ") [";
    Segment_simplex_iterator st = dt.segment_traverser_simplices_begin(s.source(), s.target());
    Segment_simplex_iterator stend = dt.segment_traverser_simplices_end();

    unsigned int inf = 0, fin = 0;
    unsigned int nb_facets = 0, nb_edges = 0, nb_vertex = 0;
    for (; st != st.end(); ++st)
    {
      std::cout << st->dimension() << " ";
      std::cout.flush();
      if (st->dimension() == 3)
      {
        if (dt.is_infinite(Cell_handle(*st))) ++inf;
        else                                  ++fin;
      }
      if (st->dimension() == 2)      ++nb_facets;
      else if (st->dimension() == 1) ++nb_edges;
      else if (st->dimension() == 0) ++nb_vertex;
    }
    std::cout << "\b]" << std::endl;

    std::cout << "\tinfinite cells : " << inf << std::endl;
    std::cout << "\tfinite cells   : " << fin << std::endl;
    std::cout << "\tfacets         : " << nb_facets << std::endl;
    std::cout << "\tedges          : " << nb_edges << std::endl;
    std::cout << "\tvertices       : " << nb_vertex << std::endl;
    std::cout << std::endl;
  }
}

int main(int argc, char* argv[])
{
  const char* fname = (argc>1) ? argv[1] : "data/blobby.xyz";
  int nb_seg = (argc > 2) ? atoi(argv[2]) : 3;

  // Reads a .xyz point set file in points.
  // As the point is the second element of the tuple (that is with index 1)
  // we use a property map that accesses the 1st element of the tuple.

  std::vector<Point_3> points;
  std::ifstream stream(fname);
  if (!stream ||
    !CGAL::IO::read_XYZ(stream, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  //bbox
  //min (-0.481293,-0.220929,-0.194076), max (0.311532,0.225525,0.198025)

  // Construct the Delaunay triangulation.
  DT dt(points.begin(), points.end());
  assert(dt.is_valid());

  CGAL::Random rng;
  for (int i = 0; i < nb_seg; ++i)
  {
    // Construct a traverser.
    Point_3 p1(rng.get_double(-0.48, 0.31),
               rng.get_double(-0.22, 0.22),
               rng.get_double(-0.19, 0.19));
    Point_3 p2(rng.get_double(-0.48, 0.31),
               rng.get_double(-0.22, 0.22),
               rng.get_double(-0.19, 0.19));

    std::cout << "Traverser " << (i + 1)
      << "\n\t(" << p1
      << ")\n\t(" << p2 << ")" << std::endl;
    Segment_simplex_iterator st = dt.segment_traverser_simplices_begin(p1, p2);
    Segment_simplex_iterator stend = dt.segment_traverser_simplices_end();

    // Count the number of finite cells traversed.
    unsigned int inf = 0, fin = 0;
    unsigned int nb_facets = 0, nb_edges = 0, nb_vertex = 0;
    for (; st != stend; ++st)
    {
      if (st->dimension() == 3)
      {
        if (Cell_handle(*st) != Cell_handle()
          && dt.is_infinite(Cell_handle(*st)))
          ++inf;
        else
          ++fin;
      }
      if (st->dimension() == 2)      ++nb_facets;
      else if (st->dimension() == 1) ++nb_edges;
      else if (st->dimension() == 0) ++nb_vertex;
    }

    std::cout << "While traversing from " << p1
              << " to " << p2 << std::endl;
    std::cout << "\tinfinite cells : " << inf << std::endl;
    std::cout << "\tfinite cells   : " << fin << std::endl;
    std::cout << "\tfacets   : " << nb_facets << std::endl;
    std::cout << "\tedges    : " << nb_edges << std::endl;
    std::cout << "\tvertices : " << nb_vertex << std::endl;
    std::cout << std::endl << std::endl;
  }

  //check degenerate cases
  // - along an edge
  test_vertex_edge_vertex(dt, 2);

  // - along a facet via edge/facet/edge
  test_edge_facet_edge(dt, 3);

  // - along a facet via edge/facet/vertex
  test_edge_facet_vertex(dt, 3);

  // - along a facet via vertex/facet/edge
  test_vertex_facet_edge(dt, 3);

  // - along 2 successive facets (vertex/facet/edge/facet/edge)
  // - along 2 successive edges (vertex/edge/vertex/edge/vertex)
  // - along a facet and an edge successively
  test_triangulation_on_a_grid();

  return 0;
}
