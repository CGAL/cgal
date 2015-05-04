//== INCLUDES =================================================================

#include <fstream>
#include "performance_2.h"
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Memory_sizer.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
typedef double Scalar;

//== CLASS DEFINITION =========================================================


class Surface_mesh_performance : public Performance_test_2
{
public:

  Surface_mesh_performance() : Performance_test_2()
  {
    points = mesh.property_map<Surface_mesh::Vertex_index,Point_3>("v:point").first;
  }


private:

  Surface_mesh mesh;
  Surface_mesh::Property_map<Surface_mesh::Vertex_index,Point_3> points;
  Surface_mesh::Property_map<Surface_mesh::Vertex_index,Vector_3> vnormals;
  Surface_mesh::Property_map<Surface_mesh::Face_index,Vector_3> fnormals;


private:

  void display_info()
  {
    std::cout << "#Darts=" << mesh.num_halfedges();
    std::cout << ", #0-cells=" << mesh.num_vertices();
    std::cout << ", #1-cells=" << mesh.num_edges();
    std::cout << ", #2-cells=" << mesh.num_faces();
    std::cout << "\t" << std::endl;
  }

  virtual bool read_mesh(const char* _filename)
  {
    CGAL::Memory_sizer ms;
    bool b = CGAL::read_off(mesh, _filename);
    std::cout << "memory consumption: " << ms.virtual_size() << "  " << ms.resident_size() << std::endl;
    return b;
  }


  virtual bool write_mesh(const char* _filename)
  {
    return CGAL::write_off(mesh, _filename);
  }


  virtual int circulator_test()
  {
    Surface_mesh::Vertex_iterator vit, vend=mesh.vertices_end();
    Surface_mesh::Face_iterator   fit, fend=mesh.faces_end();

    Surface_mesh::Face_around_target_circulator   vfit, vfend;
    Surface_mesh::Vertex_around_face_circulator   fvit, fvend;

    int counter = 0;

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
    {
      vfit = vfend = Surface_mesh::Face_around_target_circulator(mesh.halfedge(*vit),mesh);
      if (vfit) do
                {
                  ++counter;
                }
        while (++vfit != vfend);
    }

    for (fit=mesh.faces_begin(); fit!=fend; ++fit)
    {
      fvit = fvend = Surface_mesh::Vertex_around_face_circulator(mesh.halfedge(*fit),mesh);
      do
      {
        --counter;
      }
      while (++fvit != fvend);
    }

    return counter;
  }


  virtual void barycenter_test(bool draw)
  {
    Vector_3 p(0,0,0);
    Surface_mesh::Vertex_iterator vit, vend=mesh.vertices_end();

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
      p = p + (points[*vit] - CGAL::ORIGIN);

    p = p / mesh.num_vertices();

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
      points[*vit] = points[*vit] - p;

    if ( draw ) std::cout<<"Barycenter: "<<p<<std::endl;
  }


  virtual void normal_test()
  {

    vnormals = mesh.add_property_map<Surface_mesh::Vertex_index,Vector_3>("v:normal").first;
    fnormals = mesh.add_property_map<Surface_mesh::Face_index,Vector_3>("f:normal").first;
    Surface_mesh::Vertex_iterator vit, vend=mesh.vertices_end();
    Surface_mesh::Face_iterator   fit, fend=mesh.faces_end();
    Surface_mesh::Face_around_target_circulator vfit, vfend;

    for (fit=mesh.faces_begin(); fit!=fend; ++fit)
    {
      Surface_mesh::Halfedge_index  h = mesh.halfedge(*fit);
      const Point_3& p0 = points[mesh.target(h)];
      h = mesh.next(h);
      const Point_3& p1 = points[mesh.target(h)];
      h = mesh.next(h);
      const Point_3& p2 = points[mesh.target(h)];
      Vector_3 n = cross_product(p0-p1, p2-p1);
      n = n / sqrt(n.squared_length());
      fnormals[*fit] = n;
    }

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
    {
      Vector_3 n(0,0,0);
      vfit = vfend = Surface_mesh::Face_around_target_circulator(mesh.halfedge(*vit),mesh);
      if (vfit) do
                {
                  n = n + fnormals[*vfit];
                }
        while (++vfit != vfend);
      n = n / sqrt(n.squared_length());
      vnormals[*vit] = n;
    }
    mesh.remove_property_map(vnormals);
    mesh.remove_property_map(fnormals);
  }


  virtual void smoothing_test()
  {
    Surface_mesh::Vertex_iterator vit, vend=mesh.vertices_end();
    Surface_mesh::Vertex_around_target_circulator vvit, vvend;

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
    {
      if (!mesh.is_border(*vit))
      {
        Vector_3  p(0,0,0);
        Scalar c(0);
        vvit = vvend = Surface_mesh::Vertex_around_target_circulator(mesh.halfedge(*vit),mesh);
        do
        {
          p = p + (points[*vvit] - CGAL::ORIGIN);
          ++c;
        } while (++vvit != vvend);
        p = p/c;
        points[*vit] = CGAL::ORIGIN + p;
      }
    }
  }

  virtual void subdivision_test()
  {
    // reserve memory
    int nv = mesh.num_vertices();
    int ne = mesh.num_edges();
    int nf = mesh.num_faces();
    mesh.reserve(nv+nf, ne+3*nf, 3*nf);


    // iterators
    Surface_mesh::Vertex_iterator vit, vend=mesh.vertices_end();
    Surface_mesh::Face_iterator fit, fend=mesh.faces_end();
    Surface_mesh::Edge_iterator eit, eend=mesh.edges_end();


    // compute new positions of old vertices
    Surface_mesh::Property_map<Surface_mesh::Vertex_index,Point_3> new_pos = mesh.add_property_map<Surface_mesh::Vertex_index,Point_3>("v:np").first;
    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
    {
      if (!mesh.is_border(*vit))
      {
        Scalar n = mesh.degree(*vit);
        Scalar alpha = (4.0 - 2.0*cos(2.0*CGAL_PI/n)) / 9.0;
        Vector_3  p(0,0,0);
        Surface_mesh::Vertex_around_target_circulator vvit(mesh.halfedge(*vit),mesh), vvend(vvit);
        do
        {
          p = p + (points[*vvit] - CGAL::ORIGIN);
        } while (++vvit != vvend);

        p =  (1.0f-alpha)*(points[*vit] - CGAL::ORIGIN) + alpha/n*p;
        new_pos[*vit] = CGAL::ORIGIN + p;
      }
    }

    // split faces
    for (fit=mesh.faces_begin(); fit!=fend; ++fit)
    {
      Point_3  p(0,0,0);
      Scalar c(0);
      Surface_mesh::Vertex_around_face_circulator fvit(mesh.halfedge(*fit),mesh), fvend(fvit);
      do
      {
        p = p + (points[*fvit] - CGAL::ORIGIN);
        ++c;
      } while (++fvit!=fvend);

      p = CGAL::ORIGIN + (p-CGAL::ORIGIN)/c;
      Surface_mesh::Vertex_index v = mesh.target(CGAL::Euler::add_center_vertex(mesh.halfedge(*fit),mesh));
      points[v] = p;
    }

    // set new positions of old vertices
    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
      if (!mesh.is_border(*vit))
        points[*vit] = new_pos[*vit];
    mesh.remove_property_map(new_pos);

    // flip old edges
    for (eit=mesh.edges_begin(); eit!=eend; ++eit){
      if (! mesh.is_border(*eit)){
        CGAL::Euler::flip_edge(mesh.halfedge(*eit),mesh);
      }
    }
  }


  virtual void collapse_test()
  {
    // reserve memory
    int nv = mesh.num_vertices();
    int ne = mesh.num_edges();
    int nf = mesh.num_faces();
    mesh.reserve(nv+nf, ne+3*nf, 3*nf);


    // iterators
    Surface_mesh::Vertex_iterator vit, vend=mesh.vertices_end();
    Surface_mesh::Face_iterator fit, fend=mesh.faces_end();
    // split faces
    Point_3  p(0,0,0);
    for (fit=mesh.faces_begin(); fit!=fend; ++fit){
      Surface_mesh::Vertex_index v = mesh.target(CGAL::Euler::add_center_vertex(mesh.halfedge(*fit),mesh));
      points[v] = p;
    }
    int i = 0;
    // collapse new edges
    vit = vend; vend=mesh.vertices_end();
    for (; vit!=vend; ++vit){
      Surface_mesh::Halfedge_index he = mesh.halfedge(*vit);
      Surface_mesh::Vertex_index vd = mesh.source(he);
      Point_3 p = mesh.point(vd);
      Surface_mesh::Vertex_index vkept = CGAL::Euler::collapse_edge(mesh.edge(he), mesh);
      mesh.point(vkept) = p;
    }
    // remove deleted items
    //mesh.collect_garbage();
  }


  virtual void lindstrom_test(const char* _filename)
  {
    namespace SMS = CGAL::Surface_mesh_simplification ;

    mesh.clear();
    bool b = CGAL::read_off(mesh, _filename);
    SMS::Count_ratio_stop_predicate<Surface_mesh> stop(0.1);
    int r = SMS::edge_collapse(mesh, stop);
  }


};
//=============================================================================
