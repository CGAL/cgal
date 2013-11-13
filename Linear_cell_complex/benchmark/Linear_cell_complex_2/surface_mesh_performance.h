//== INCLUDES =================================================================


#include "performance_2.h"
#include "surface_mesh/Surface_mesh.h"


//== CLASS DEFINITION =========================================================


class Surface_mesh_performance : public Performance_test_2
{
public:

  Surface_mesh_performance() : Performance_test_2()
  {
    points   = mesh.vertex_property<Point>("v:point");
    vnormals = mesh.vertex_property<Point>("v:normal");
    fnormals = mesh.face_property<Point>("f:normal");
  }


private:

  Surface_mesh mesh;
  Surface_mesh::Vertex_property<Point> points;
  Surface_mesh::Vertex_property<Point> vnormals;
  Surface_mesh::Face_property<Point>   fnormals;


private:

  void display_info()
  {
    std::cout << "#Darts=" << mesh.n_halfedges();
    std::cout << ", #0-cells=" << mesh.n_vertices();
    std::cout << ", #1-cells=" << mesh.n_edges();
    std::cout << ", #2-cells=" << mesh.n_faces();
    std::cout << "\t" << std::endl;
  }

  virtual bool read_mesh(const char* _filename)
  {
    return mesh.read(_filename);
  }


  virtual bool write_mesh(const char* _filename)
  {
    return mesh.write(_filename);
  }


  virtual int circulator_test()
  {
    Surface_mesh::Vertex_iterator vit, vend=mesh.vertices_end();
    Surface_mesh::Face_iterator   fit, fend=mesh.faces_end();

    Surface_mesh::Face_around_vertex_circulator   vfit, vfend;
    Surface_mesh::Vertex_around_face_circulator   fvit, fvend;

    int counter = 0;

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
    {
      vfit = vfend = mesh.faces(*vit);
      if (vfit) do
                {
                  ++counter;
                }
        while (++vfit != vfend);
    }

    for (fit=mesh.faces_begin(); fit!=fend; ++fit)
    {
      fvit = fvend = mesh.vertices(*fit);
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
    Point p(0,0,0);
    Surface_mesh::Vertex_iterator vit, vend=mesh.vertices_end();

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
      p += points[*vit];

    p /= mesh.n_vertices();

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
      points[*vit] -= p;

    if ( draw ) std::cout<<"Barycenter: "<<p<<std::endl;
  }


  virtual void normal_test()
  {
    Surface_mesh::Vertex_iterator vit, vend=mesh.vertices_end();
    Surface_mesh::Face_iterator   fit, fend=mesh.faces_end();
    Surface_mesh::Face_around_vertex_circulator vfit, vfend;

    for (fit=mesh.faces_begin(); fit!=fend; ++fit)
    {
      Surface_mesh::Halfedge h = mesh.halfedge(*fit);
      Point p0 = points[mesh.to_vertex(h)];
      h = mesh.next_halfedge(h);
      Point p1 = points[mesh.to_vertex(h)];
      p1 -= p0;
      h = mesh.next_halfedge(h);
      Point p2 = points[mesh.to_vertex(h)];
      p2 -= p0;
      fnormals[*fit] = cross(p1, p2).normalize();
    }

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
    {
      Point n(0,0,0);
      vfit = vfend = mesh.faces(*vit);
      if (vfit) do
                {
                  n += fnormals[*vfit];
                }
        while (++vfit != vfend);
      vnormals[*vit] = n.normalize();
    }
  }


  virtual void smoothing_test()
  {
    Surface_mesh::Vertex_iterator vit, vend=mesh.vertices_end();
    Surface_mesh::Vertex_around_vertex_circulator vvit, vvend;

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
    {
      if (!mesh.is_boundary(*vit))
      {
        Point  p(0,0,0);
        Scalar c(0);
        vvit = vvend = mesh.vertices(*vit);
        do
        {
          p += points[*vvit];
          ++c;
        }
        while (++vvit != vvend);
        p /= c;
        points[*vit] = p;
      }
    }
  }


  virtual void subdivision_test()
  {
    // reserve memory
    int nv = mesh.n_vertices();
    int ne = mesh.n_edges();
    int nf = mesh.n_faces();
    mesh.reserve(nv+nf, ne+3*nf, 3*nf);


    // iterators
    Surface_mesh::Vertex_iterator vit, vend=mesh.vertices_end();
    Surface_mesh::Face_iterator fit, fend=mesh.faces_end();
    Surface_mesh::Edge_iterator eit, eend=mesh.edges_end();


    // compute new positions of old vertices
    Surface_mesh::Vertex_property<Point> new_pos = mesh.add_vertex_property<Point>("v:np");
    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
    {
      if (!mesh.is_boundary(*vit))
      {
        Scalar n = mesh.valence(*vit);
        Scalar alpha = (4.0 - 2.0*cos(2.0*M_PI/n)) / 9.0;
        Point  p(0,0,0);
        Surface_mesh::Vertex_around_vertex_circulator vvit=mesh.vertices(*vit), vvend=vvit;
        do
        {
          p += points[*vvit];
        }
        while (++vvit != vvend);
        p = (1.0f-alpha)*points[*vit] + alpha/n*p;
        new_pos[*vit] = p;
      }
    }


    // split faces
    for (fit=mesh.faces_begin(); fit!=fend; ++fit)
    {
      Point  p(0,0,0);
      Scalar c(0);
      Surface_mesh::Vertex_around_face_circulator fvit = mesh.vertices(*fit), fvend=fvit;
      do
      {
        p += points[*fvit];
        ++c;
      }
      while (++fvit!=fvend);
      p /= c;

      mesh.split(*fit, p);
    }


    // set new positions of old vertices
    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
      if (!mesh.is_boundary(*vit))
        points[*vit] = new_pos[*vit];
    mesh.remove_vertex_property(new_pos);


    // flip old edges
    for (eit=mesh.edges_begin(); eit!=eend; ++eit)
      if (mesh.is_flip_ok(*eit))
        mesh.flip(*eit);
  }


  virtual void collapse_test()
  {
    // reserve memory
    int nv = mesh.n_vertices();
    int ne = mesh.n_edges();
    int nf = mesh.n_faces();
    mesh.reserve(nv+nf, ne+3*nf, 3*nf);


    // iterators
    Surface_mesh::Vertex_iterator vit, vend=mesh.vertices_end();
    Surface_mesh::Face_iterator fit, fend=mesh.faces_end();

    // split faces
    Point  p(0,0,0);
    for (fit=mesh.faces_begin(); fit!=fend; ++fit)
      mesh.split(*fit, p);

    // collapse new edges
    vit = vend; vend=mesh.vertices_end();
    for (; vit!=vend; ++vit)
      mesh.collapse(mesh.halfedge(*vit));

    // remove deleted items
    mesh.garbage_collection();
  }

};
//=============================================================================
