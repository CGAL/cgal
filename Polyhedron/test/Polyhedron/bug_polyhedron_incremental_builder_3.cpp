#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <sstream>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef Polyhedron::HalfedgeDS HDS;
typedef CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
typedef CGAL::HalfedgeDS_decorator<HDS> Decorator;
typedef HDS::Vertex::Base VBase;
struct Appender : public CGAL::Modifier_base<HDS>
{
  void operator()(HDS& hds)
  {
    Builder B(hds);
    B.begin_surface(1,1,4,Builder::ABSOLUTE_INDEXING);
    B.add_vertex( K::Point_3(2,2,0) );
    B.begin_facet();
    B.add_vertex_to_facet(3);
    B.add_vertex_to_facet(4);
    B.add_vertex_to_facet(6);
    B.end_facet();
    B.end_surface();
  }
};

/*!
We start with this configuration:
  0---3---4
  |\  |\  |
  | \ | \ |
  |  \|  \|
  1---2---5

Then we update a the halfedge pointer of vertices
3 and 4 so that the halfedges are on the boundary.
We append a new facet to get the following:
          6
         /|
        / |
       /  |
  0---3---4
  |\  |\  |
  | \ | \ |
  |  \|  \|
  1---2---5
*/

int main()
{
  std::stringstream ss;
  ss << "\
OFF\n6 4 0\n\
0 1 0\n\
0 0 0\n\
1 0 0\n\
1 1 0\n\
2 1 0\n\
2 0 0\n\
3 0 1 2\n\
3 0 2 3\n\
3 2 5 3\n\
3 5 4 3\n";

  Polyhedron P;
  ss >> P;

  assert( P.size_of_vertices() == 6);
  assert( P.size_of_facets() == 4);
  assert( P.is_valid() );

  //consider vertex 3 and set its halfedge to be on the border
  Polyhedron::Vertex_iterator vit=P.vertices_begin();
  std::advance(vit, 3);
  assert( vit->point() == K::Point_3(1, 1, 0) );
  Polyhedron::Halfedge_handle h=vit->halfedge();
  while ( !h->is_border() )
    h=h->next()->opposite();
  vit->VBase::set_halfedge(h);
  assert( vit->halfedge()->vertex() == vit );

  //consider vertex 4 and set its halfedge to be on the border
  ++vit;
  assert( vit->point() == K::Point_3(2, 1, 0) );
  h=vit->halfedge();
  while ( !h->is_border() )
    h=h->next()->opposite();
  vit->VBase::set_halfedge(h);
  assert( vit->halfedge()->vertex() == vit );

  //try to append a facet
  Appender modifier;
  P.delegate(modifier);

  assert( P.size_of_vertices() == 7);
  assert( P.size_of_facets() == 5);
  assert( P.is_valid() );
}
