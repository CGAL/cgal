#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef CGAL::Surface_mesh_simplification::Edge_profile<Mesh> Profile;

void naive_all_triangles(Mesh::Halfedge_index h, Mesh& m, std::set<Mesh::Face_index>& triangles)
{
  for(Mesh::Halfedge_index hh : CGAL::halfedges_around_source(h, m))
  {
    if(!is_border(hh, m))
      triangles.insert(face(hh,m));
  }
  for(Mesh::Halfedge_index hh : CGAL::halfedges_around_target(h, m))
  {
    if(!is_border(hh, m))
      triangles.insert(face(hh,m));
  }
}

void naive_link_vertices(Mesh::Halfedge_index h, Mesh& m,
                         const std::set<Mesh::Face_index>& triangles,
                         std::set<Mesh::Vertex_index>& link_vertices)
{
  for(Mesh::Face_index f : triangles)
  {
    for(Mesh::Halfedge_index h : CGAL::halfedges_around_face(halfedge(f, m), m))
    {
      link_vertices.insert(target(h, m));
    }
  }
  link_vertices.erase(source(h, m));
  link_vertices.erase(target(h, m));
}

struct A{};

boost::tuple<Mesh::Vertex_index, Mesh::Vertex_index, Mesh::Vertex_index>
make_canonical_tuple(Mesh::Vertex_index v1, Mesh::Vertex_index v2, Mesh::Vertex_index v3)
{
  Mesh::Vertex_index vs[3]={v1, v2, v3};
  std::sort(&vs[0], &vs[0]+3);
  return boost::make_tuple(vs[0], vs[1], vs[2]);
}

boost::tuple<Mesh::Vertex_index, Mesh::Vertex_index, Mesh::Vertex_index>
make_canonical_tuple(const Profile::Triangle& t)
{
  return make_canonical_tuple(t.v0, t.v1, t.v2);
}

boost::tuple<Mesh::Vertex_index, Mesh::Vertex_index, Mesh::Vertex_index>
make_canonical_tuple(Mesh::Face_index f, Mesh& m)
{
  Mesh::Halfedge_index h=halfedge(f, m);
  return make_canonical_tuple(source(h,m), target(h,m), target(next(h, m), m));
}

void test(const char* fname)
{
  Mesh m;
  std::ifstream input(fname);
  assert(!input.fail());
  input >> m;
  assert(num_vertices(m)!=0);
  A a;

  for(Mesh::Halfedge_index h : halfedges(m))
  {
    std::set<Mesh::Face_index> triangles;
    naive_all_triangles(h, m, triangles);

    Profile profile(h, m, K(), a, get(boost::vertex_point, m), a, true);

    if(CGAL::Euler::does_satisfy_link_condition(edge(h, m), m))
    {
      std::set<Mesh::Vertex_index> link_vertices;
      naive_link_vertices(h, m, triangles, link_vertices);
      assert(link_vertices.size()==profile.link().size());
      assert(std::set<Mesh::Vertex_index>(profile.link().begin(),
                                           profile.link().end()).size()
              == link_vertices.size());

      for(const Mesh::Vertex_index& v : profile.link())
      {
        assert(link_vertices.count(v) == 1);
      }
    }

    assert(triangles.size() == profile.triangles().size());
    std::set<boost::tuple<Mesh::Vertex_index, Mesh::Vertex_index, Mesh::Vertex_index> > triple_set;
    for(const Profile::Triangle& t : profile.triangles())
    {
      triple_set.insert(make_canonical_tuple(t));
    }
    for(Mesh::Face_index f : triangles)
    {
      assert(triple_set.count(make_canonical_tuple(f, m)) == 1);
    }
  }
}

int main(int argc, char** argv)
{
  for(int i=1; i<argc; ++i)
  {
    std::cout << "Testing " << argv[i] << "\n";
    test(argv[i]);
  }
  if(argc==1)
  {
    std::cout << "No file provided, nothing tested\n";
  }

  return 0;
}

