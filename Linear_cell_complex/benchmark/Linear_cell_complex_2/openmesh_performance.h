//== INCLUDES =================================================================


#include "performance_2.h"
#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"


//== CLASS DEFINITION =========================================================


class OpenMesh_performance : public Performance_test_2
{
public:

  OpenMesh_performance() : Performance_test_2()
  {
    mesh.request_vertex_status();
    mesh.request_edge_status();
    mesh.request_face_status();
    mesh.request_vertex_normals();
    mesh.request_face_normals();
  }


private:

  struct MyTraits : public OpenMesh::DefaultTraits
  {
    typedef OpenMesh::Vec3d Point;
    typedef OpenMesh::Vec3d Normal;
  };

  typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> Mesh;
  Mesh mesh;


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
    return OpenMesh::IO::read_mesh(mesh, _filename);
  }


  virtual bool write_mesh(const char* _filename)
  {
    return OpenMesh::IO::write_mesh(mesh, _filename);
  }

  virtual int circulator_test()
  {
    Mesh::VIter vit, vend=mesh.vertices_end();
    Mesh::FIter fit, fend=mesh.faces_end();
    int counter = 0;

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
      for (Mesh::VFIter vfit=mesh.vf_iter(vit); vfit; ++vfit)
        ++counter;

    for (fit=mesh.faces_begin(); fit!=fend; ++fit)
      for (Mesh::FVIter fvit=mesh.fv_iter(fit); fvit; ++fvit)
        --counter;

    return counter;
  }

  virtual void barycenter_test(bool draw)
  {
    Mesh::VIter vit, vend=mesh.vertices_end();
    Mesh::Point p(0,0,0);

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
      p += mesh.point(vit);

    p /= mesh.n_vertices();

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
      mesh.point(vit) -= p;

    if ( draw ) std::cout<<"Barycenter: "<<p<<std::endl;
  }

  virtual void normal_test()
  {
    Mesh::VIter vit, vend=mesh.vertices_end();
    Mesh::FIter fit, fend=mesh.faces_end();
    Mesh::VFIter vfit;

    for (fit=mesh.faces_begin(); fit!=fend; ++fit)
    {
      Mesh::HalfedgeHandle h = mesh.halfedge_handle(fit);
      Mesh::Point p0 = mesh.point(mesh.to_vertex_handle(h));
      h = mesh.next_halfedge_handle(h);
      Mesh::Point p1 = mesh.point(mesh.to_vertex_handle(h));
      h = mesh.next_halfedge_handle(h);
      Mesh::Point p2 = mesh.point(mesh.to_vertex_handle(h));
      mesh.set_normal(fit, ((p2-=p1)%(p0-=p1)).normalize());
    }

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
    {
      Mesh::Point n(0,0,0);
      for (vfit=mesh.vf_iter(vit); vfit; ++vfit)
        n += mesh.normal(vfit);
      mesh.set_normal(vit, n.normalize());
    }
  }

  virtual void smoothing_test()
  {
    Mesh::VIter vit, vend=mesh.vertices_end();

    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
    {
      if (!mesh.is_boundary(vit))
      {
        Mesh::Point  p(0,0,0);
        Mesh::Scalar c(0);
        for (Mesh::VVIter vvit = mesh.vv_iter(vit); vvit; ++vvit)
        {
          p += mesh.point(vvit);
          ++c;
        }
        p /= c;
        mesh.point(vit) = p;
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
    Mesh::VIter vit, vend=mesh.vertices_end();
    Mesh::FIter fit, fend=mesh.faces_end();
    Mesh::EIter eit, eend=mesh.edges_end();

    // compute new positions of old vertices
    OpenMesh::VPropHandleT<Mesh::Point> new_pos;
    mesh.add_property(new_pos);
    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
    {
      if (!mesh.is_boundary(vit))
      {
        Mesh::Scalar n = mesh.valence(vit);
        Mesh::Scalar alpha = (4.0 - 2.0*cos(2.0*M_PI/n)) / 9.0;
        Mesh::Point  p(0,0,0);
        for (Mesh::VVIter vvit = mesh.vv_iter(vit); vvit; ++vvit)
          p += mesh.point(vvit);
        p = (1.0f-alpha)*mesh.point(vit) + alpha/n*p;
        mesh.property(new_pos, vit) = p;
      }
    }

    // split faces
    for (fit=mesh.faces_begin(); fit!=fend; ++fit)
    {
      Mesh::Point  p(0,0,0);
      Mesh::Scalar c(0);
      for (Mesh::FVIter fvit=mesh.fv_iter(fit); fvit; ++fvit)
      {
        p += mesh.point(fvit);
        ++c;
      }
      p /= c;

      mesh.split(fit, p);
    }

    // set new positions of old vertices
    for (vit=mesh.vertices_begin(); vit!=vend; ++vit)
      if (!mesh.is_boundary(vit))
        mesh.point(vit) = mesh.property(new_pos, vit);
    mesh.remove_property(new_pos);


    // flip old edges
    for (eit=mesh.edges_begin(); eit!=eend; ++eit)
      if (mesh.is_flip_ok(eit))
        mesh.flip(eit);
  }

  virtual void collapse_test()
  {
    // reserve memory
    int nv = mesh.n_vertices();
    int ne = mesh.n_edges();
    int nf = mesh.n_faces();
    mesh.reserve(nv+nf, ne+3*nf, 3*nf);

    // iterators
    Mesh::VIter vit, vend=mesh.vertices_end();
    Mesh::FIter fit, fend=mesh.faces_end();

    // split faces
    Mesh::Point p(0,0,0);
    for (fit=mesh.faces_begin(); fit!=fend; ++fit)
      mesh.split(fit, p);

    // collapse new edges
    vit = vend; vend=mesh.vertices_end();
    for (; vit!=vend; ++vit)
      mesh.collapse(mesh.halfedge_handle(vit));

    // remove deleted items
    mesh.garbage_collection();
  }

  virtual void remesh_test()
  {

  }
};
//=============================================================================
