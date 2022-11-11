// Use the default kernels and the default conversions:

// *********************
// Regular triangulation
// *********************
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cassert>
#include <CGAL/Regular_triangulation_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_3<K>       Regular;

#include <CGAL/Compute_anchor_3.h>

#include <fstream>

#include <CGAL/Triangulation_simplex_3.h>
typedef K::RT                                  Weight;
typedef K::Point_3                             Point;
typedef K::Weighted_point_3                    Weighted_point;
typedef Regular::Vertex_handle                 Vertex_handle;
typedef Regular::Edge                          Edge;
typedef Regular::Facet                         Facet;
typedef Regular::Cell_handle                   Cell_handle;
typedef CGAL::Triangulation_simplex_3<Regular> Simplex;

void test_anchor_del()
{
  Regular reg;
  CGAL::Compute_anchor_3<Regular> compute_anchor_obj(reg);
  Vertex_handle vh[4];
  Cell_handle ch;
  Simplex s;
  int index1, index2;

  vh[0] = reg.insert(Weighted_point(Point(-1,0,0), 1.1));
  vh[1] = reg.insert(Weighted_point(Point( 1,0,0), 1.1));
  vh[2] = reg.insert(Weighted_point(Point( 0,1,0), 1));
  vh[3] = reg.insert(Weighted_point(Point( 0,0,1), 1));
  ch = reg.finite_cells_begin();

  // Regular case:
  for (int i=0; i<4; i++) {
    assert(Simplex(vh[i]) == compute_anchor_obj.anchor_del(vh[i]));
    assert(Simplex(vh[i]) == compute_anchor_obj.anchor_vor(vh[i]));
    for (int j=i+1; j<4; j++) {
      assert(Simplex(Edge(ch,i,j)) ==
             compute_anchor_obj.anchor_del(Edge(ch,i,j)));
      assert(Simplex(Edge(ch,i,j)) ==
             compute_anchor_obj.anchor_vor(Edge(ch,i,j)));
    }
    assert(Simplex(Facet(ch,i)) ==
           compute_anchor_obj.anchor_del(Facet(ch,i)));
    assert(Simplex(Facet(ch,i)) ==
           compute_anchor_obj.anchor_vor(Facet(ch,i)));
    assert(Simplex(ch) == compute_anchor_obj.anchor_del(ch));
    assert(Simplex(ch) == compute_anchor_obj.anchor_vor(ch));
  }
  reg.clear();

  vh[0] = reg.insert(Weighted_point(Point(-1,0,0), 1));
  vh[1] = reg.insert(Weighted_point(Point( 1,0,0), 1));
  vh[2] = reg.insert(Weighted_point(Point( 0,1,0), 1));
  vh[3] = reg.insert(Weighted_point(Point( 0,0,1), 1));
  ch = reg.finite_cells_begin();

  // Voronoi vertex on a Delaunay edge (degenerate):
  for (int i=0; i<4; i++) {
    assert(Simplex(vh[i]) == compute_anchor_obj.anchor_del(vh[i]));
    assert(Simplex(vh[i]) == compute_anchor_obj.anchor_vor(vh[i]));
    assert(Simplex(Facet(ch,i)) ==
           compute_anchor_obj.anchor_del(Facet(ch,i)));
    assert(Simplex(Facet(ch,i)) ==
           compute_anchor_obj.anchor_vor(Facet(ch,i)));
    for (int j=i+1; j<4; j++) {
      assert(Simplex(Edge(ch,i,j)) ==
             compute_anchor_obj.anchor_del(Edge(ch,i,j)));
      assert(Simplex(Edge(ch,i,j)) ==
             compute_anchor_obj.anchor_vor(Edge(ch,i,j)));
    }
    assert(Simplex(ch) == compute_anchor_obj.anchor_del(ch));
    assert(Simplex(ch) == compute_anchor_obj.anchor_vor(ch));
  }
  reg.clear();

  vh[0] = reg.insert(Weighted_point(Point(0,0,0), .99));
  vh[1] = reg.insert(Weighted_point(Point(1,0,0), 2));
  vh[2] = reg.insert(Weighted_point(Point(0,1,0), 2));
  vh[3] = reg.insert(Weighted_point(Point(0,0,1), 2));
  ch = reg.finite_cells_begin();
  index1 = ch->index(vh[0]);

  // Anchor of a Delaunay edge/facet/cell on a vertex
  assert(Simplex(vh[0]) ==
      compute_anchor_obj.anchor_del(Edge(ch,index1,(index1+1)&3)));
  assert(Simplex(vh[0]) ==
      compute_anchor_obj.anchor_del(Edge(ch,index1,(index1+2)&3)));
  assert(Simplex(vh[0]) ==
      compute_anchor_obj.anchor_del(Edge(ch,index1,(index1+3)&3)));
  assert(Simplex(vh[0]) ==
      compute_anchor_obj.anchor_del(Facet(ch,(index1+1)&3)));
  assert(Simplex(vh[0]) ==
      compute_anchor_obj.anchor_del(Facet(ch,(index1+2)&3)));
  assert(Simplex(vh[0]) ==
      compute_anchor_obj.anchor_del(Facet(ch,(index1+3)&3)));
  assert(Simplex(vh[0]) == compute_anchor_obj.anchor_del(ch));
  reg.clear();

  vh[0] = reg.insert(Weighted_point(Point(-1,0,0), 1));
  vh[1] = reg.insert(Weighted_point(Point( 1,0,0), 1));
  vh[2] = reg.insert(Weighted_point(Point( 0,1,0), 1.1));
  vh[3] = reg.insert(Weighted_point(Point( 0,0,1), 1.1));
  ch = reg.finite_cells_begin();
  index1 = ch->index(vh[0]);
  index2 = ch->index(vh[1]);

  // Anchor of a Delaunay facet/cell on an edge
  assert(Simplex(Edge(ch,index1,index2)) ==
         compute_anchor_obj.anchor_del(Facet(ch,ch->index(vh[2]))));
  assert(Simplex(Edge(ch,index1,index2)) ==
         compute_anchor_obj.anchor_del(Facet(ch,ch->index(vh[3]))));
  assert(Simplex(Edge(ch,index1,index2)) ==
         compute_anchor_obj.anchor_del(ch));
  reg.clear();

  vh[0] = reg.insert(Weighted_point(Point(0,0,0), 1));
  vh[1] = reg.insert(Weighted_point(Point(1,0,0), 1));
  vh[2] = reg.insert(Weighted_point(Point(0,1,0), 1));
  vh[3] = reg.insert(Weighted_point(Point(0,0,1), 1));
  ch = reg.finite_cells_begin();
  index1 = ch->index(vh[0]);

  // Anchor of a Delaunay cell on an facet
  assert(Simplex(Facet(ch,index1)) == compute_anchor_obj.anchor_del(ch));
  reg.clear();
}

int main(int, char **)
{
  test_anchor_del();

//   Regular regular;
//   double shrink = .75;

//   {
//     std::fstream in(argv[1], std::ios::in);
//     if (!in.is_open()) {
//       std::cerr << "Could not open data file" << std::endl;
//       exit(1);
//     }
//     double xExtr[2], yExtr[2], zExtr[2];

//     double x,y,z,w;
//     in >> x >> y >> z >> w;
//     regular.insert(Weighted_point(Point(x,y,z),w));
//     if (w<0) w=0.01;
//     for (int i=0; i<2; i++) {
//       xExtr[i] = x + (-1+2*i)*std::sqrt(w);
//       yExtr[i] = y + (-1+2*i)*std::sqrt(w);
//       zExtr[i] = z + (-1+2*i)*std::sqrt(w);
//     }

//     while (in >> x >> y >> z >> w) {
//       //std::cerr << x << " " << y << " " << z << " " << w << std::endl;
//       regular.insert(Weighted_point(Point(x,y,z),w));

//       if (w<0) w=0.01;
//       xExtr[0] = std::min(xExtr[0], x - std::sqrt(w));
//       yExtr[0] = std::min(yExtr[0], y - std::sqrt(w));
//       zExtr[0] = std::min(zExtr[0], z - std::sqrt(w));
//       xExtr[1] = std::max(xExtr[1], x + std::sqrt(w));
//       yExtr[1] = std::max(yExtr[1], y + std::sqrt(w));
//       zExtr[1] = std::max(zExtr[1], z + std::sqrt(w));
//     }

//     // At least 1*max
//     double boxSize = 1 + 2.05/(shrink * (1-shrink))*
//       std::max((xExtr[1]-xExtr[0]),
//         std::max((yExtr[1]-yExtr[0]),(zExtr[1]-zExtr[0])));

//     regular.insert(Weighted_point(Point(xExtr[0]-boxSize,0,0),-1));
//     regular.insert(Weighted_point(Point(xExtr[1]+boxSize,0,0),-1));
//     regular.insert(Weighted_point(Point(0,yExtr[0]-boxSize,0),-1));
//     regular.insert(Weighted_point(Point(0,yExtr[1]+boxSize,0),-1));
//     regular.insert(Weighted_point(Point(0,0,zExtr[0]-boxSize),-1));
//     regular.insert(Weighted_point(Point(0,0,zExtr[1]+boxSize),-1));

//     assert( regular.is_valid() );
//     assert( regular.dimension() == 3 );
//   }

//   // Regular triangulation read and bounding box created

//   Simplicial simplicial;
//   Mixed_complex_builder mixed_builder(regular, simplicial, shrink);

//   assert(simplicial.is_valid(true));

  //         // NGHK: function
  //         Mesh mesh;
  //         Marching_tetrahedra marching_tetrahedra(simplicial, mesh);
  //         marching_tetrahedra.build();

  //         {
  //                 std::ofstream out("out/coarse.off");
  //                 out << mesh;
  //         }

  //         // subdivision engine
  //         Sqrt3Method subdivider(simplicial, mesh);
  //         subdivider.subdivide(2);

  //         {
  //                 std::ofstream out("out/sqrt.off");
  //                 out << mesh;
  //         }
}
