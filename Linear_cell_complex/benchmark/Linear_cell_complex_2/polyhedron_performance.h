//== INCLUDES =================================================================

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/HalfedgeDS_list.h>
#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_halfedge_base.h>
#include <CGAL/HalfedgeDS_face_base.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <fstream>
#include "performance_2.h"


//== CLASS DEFINITION =========================================================


//*****************************************************************************
// control whether Polyhedron stores vertex and face normals
#define HAS_NORMALS 1
//*****************************************************************************


typedef CGAL::Simple_cartesian<double>  CGALKernel;
typedef CGALKernel::Point_3             Point_3;
typedef CGAL::Vector_3<CGALKernel>      Vector_3;


template <class Refs>
struct MyVertex : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point_3>
{
  typedef CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point_3> Base;
  MyVertex() : Base() {}
  MyVertex(const Point_3& p) : Base(p) {}
  Vector_3 normal;
};


template <class Refs>
struct MyFace : public CGAL::HalfedgeDS_face_base<Refs>
{
  Vector_3 normal;
};


class MyItems
{
public:

  template <class Refs, class Traits>
  struct Vertex_wrapper
  {
    typedef typename Traits::Point_3 Point;
#if HAS_NORMALS
    typedef MyVertex<Refs> Vertex;
#else
    typedef CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> Vertex;
#endif
  };

  template <class Refs, class Traits>
  struct Halfedge_wrapper
  {
    typedef CGAL::HalfedgeDS_halfedge_base<Refs> Halfedge;
  };

  template <class Refs, class Traits>
  struct Face_wrapper
  {
#if HAS_NORMALS
    typedef MyFace<Refs> Face;
#else
    typedef CGAL::HalfedgeDS_face_base<Refs> Face;
#endif
  };
};

typedef CGAL::Polyhedron_3<CGALKernel, MyItems, CGAL::HalfedgeDS_list> Polyhedron;




//== CLASS DEFINITION =========================================================


class Polyhedron_performance : public Performance_test_2
{

private:

  Polyhedron P;

private:
  virtual void display_info()
  {
    std::cout<<"#darts: "<<P.size_of_halfedges()
             <<"  #0-cells: "<<P.size_of_vertices()
             <<"  #2-cells: "<<P.size_of_facets()<<std::endl;
  }

  virtual bool read_mesh(const char* _filename)
  {
    std::ifstream ifs(_filename);
    ifs >> P;
    return true;
  }


  virtual bool write_mesh(const char* _filename)
  {
    std::ofstream ofs(_filename);
    ofs << P;
    return true;
  }


  virtual int circulator_test()
  {
    Polyhedron::Vertex_iterator vit, vend=P.vertices_end();
    Polyhedron::Face_iterator   fit, fend=P.facets_end();

    Polyhedron::Halfedge_around_vertex_circulator vhit, vhend;
    Polyhedron::Halfedge_around_facet_circulator  fhit, fhend;

    int counter = 0;

    for (vit = P.vertices_begin(); vit != vend; ++vit)
    {
      vhit = vhend = vit->vertex_begin();
      do
      {
        if (!vhit->is_border())
          ++counter;
      }
      while (++vhit != vhend);
    }

    for (fit = P.facets_begin(); fit != fend; ++fit)
    {
      fhit = fhend = fit->facet_begin();
      do
      {
        --counter;
      }
      while (++fhit != fhend);
    }

    return counter;
  }


  virtual void barycenter_test(bool draw)
  {
    Polyhedron::Vertex_iterator vit, vend=P.vertices_end();
    Vector_3 v(CGAL::NULL_VECTOR);
    for (vit = P.vertices_begin(); vit != vend; ++vit)
    {
      v = v + (vit->point() - CGAL::ORIGIN);
    }
    v = v / P.size_of_vertices();
    for (vit = P.vertices_begin(); vit != vend; ++vit)
    {
      vit->point() = vit->point() - v;
    }
    if (draw) std::cout<<"Barycenter="<<v<<std::endl;
  }


  virtual void normal_test()
  {
#if HAS_NORMALS
    Polyhedron::Vertex_iterator vit, vend=P.vertices_end();
    Polyhedron::Face_iterator   fit, fend=P.facets_end();
    Polyhedron::Halfedge_around_vertex_circulator vhit, vhend;

    for (fit = P.facets_begin(); fit != fend; ++fit)
    {
      Polyhedron::Halfedge_descriptor h = fit->halfedge();
      Point_3& p0 = h->vertex()->point();
      h = h->next();
      Point_3& p1 = h->vertex()->point();
      h = h->next();
      Point_3& p2 = h->vertex()->point();
      Vector_3 n = cross_product(p0-p1, p2-p1);
      n = n / sqrt(n.squared_length());
      fit->normal = n;
    }

    for (vit=P.vertices_begin(); vit!=vend; ++vit)
    {
      Vector_3 n(0,0,0);
      vhit = vhend = vit->vertex_begin();
      do
      {
        if (!vhit->is_border())
          n = n + vhit->face()->normal;
      }
      while (++vhit != vhend);
      n = n / sqrt(n.squared_length());
      vit->normal = n;
    }
#endif
  }


  virtual void smoothing_test()
  {
    Polyhedron::Vertex_iterator vit, vend=P.vertices_end();

    for (vit = P.vertices_begin(); vit != vend; ++vit)
    {
      bool vertex_is_border = false;
      Vector_3 v(0,0,0);
      float c(0);

      Polyhedron::Halfedge_around_vertex_circulator hc = vit->vertex_begin();
      do
      {
        if (hc->is_border() || hc->opposite()->is_border())
        {
          vertex_is_border = true;
          break;
        }

        v = v + (hc->opposite()->vertex()->point() - CGAL::ORIGIN);
        ++c;
      }
      while (++hc != vit->vertex_begin());

      if (!vertex_is_border)
        vit->point() = CGAL::ORIGIN + (v / c);
    }
  }


  virtual void subdivision_test()
  {
    int nv = P.size_of_vertices();
    int nh = P.size_of_halfedges();
    int nf = P.size_of_facets();
    P.reserve(nv+nf, nh+6*nf, 3*nf);


    // iterators
    Polyhedron::Vertex_iterator vit, vend = P.vertices_end();
    Polyhedron::Face_iterator   fit, fend = P.facets_end();
    Polyhedron::Edge_iterator   eit, eend = P.edges_end();


    // compute new positions of old vertices
    int i;
    std::vector<Point_3> new_pos(nv);
    for (vit=P.vertices_begin(), i=0; vit!=vend; ++vit, ++i)
    {
      bool      is_border = false;
      Vector_3  v(0,0,0);
      float     n = CGAL::circulator_size(vit->vertex_begin());
      float     alpha = (4.0 - 2.0*cos(2.0*M_PI/n)) / 9.0;

      Polyhedron::Halfedge_around_vertex_circulator hc = vit->vertex_begin();
      do
      {
        if (hc->is_border() || hc->opposite()->is_border())
        {
          is_border = true;
          break;
        }

        v = v + (hc->opposite()->vertex()->point() - CGAL::ORIGIN);
      }
      while (++hc != vit->vertex_begin());

      v = (1.0f-alpha)*(vit->point() - CGAL::ORIGIN) + alpha/n*v;

      new_pos[i] = (is_border ? vit->point() : CGAL::ORIGIN + v);
    }


    // adjust end iterators
    --vend; --eend; --fend;


    // split faces (a for loop does not work for list-kernel!)
    fit = P.facets_begin();
    do
    {
      Vector_3  v(CGAL::NULL_VECTOR);
      float     c(0);

      Polyhedron::Halfedge_around_facet_circulator hc = fit->facet_begin();
      Polyhedron::Halfedge_around_facet_circulator hc_end(hc);
      do
      {
        v = v + (hc->vertex()->point() - CGAL::ORIGIN);
        ++c;
      }
      while (++hc != hc_end);
      v = v / c;

      Polyhedron::Halfedge_descriptor h = P.create_center_vertex(fit->halfedge());
      h->vertex()->point() = CGAL::ORIGIN + v;
    }
    while (fit++ != fend);


    // adjust end iterators
    ++vend; ++eend; ++fend;


    // set new positions of old vertices
    for (vit=P.vertices_begin(), i=0; vit!=vend; ++vit, ++i)
      vit->point() = new_pos[i];


    // flip old edges
    for (eit = P.edges_begin(); eit!=eend; ++eit)
    {
      // careful: eit->is_border() does not work
      Polyhedron::Halfedge_descriptor h = eit;
      if (!(h->is_border() || h->opposite()->is_border()))
        P.flip_edge(h);
    }
  }



  virtual void collapse_test()
  {

    // reserve memory
    int nv = P.size_of_vertices();
    int nh = P.size_of_halfedges();
    int nf = P.size_of_facets();
    P.reserve(nv+nf, nh+6*nf, 3*nf);

    // iterators
    Polyhedron::Vertex_iterator vit, vend = P.vertices_end();
    Polyhedron::Face_iterator   fit, fend = P.facets_end();


    // adjust end iterators
    --vend; --fend;


    // split faces (a for loop does not work for list-kernel!)
    Point_3  p(0,0,0);
    fit = P.facets_begin();
    do
    {
      Polyhedron::Halfedge_descriptor h = P.create_center_vertex(fit->halfedge());
      h->vertex()->point() = p;
    }
    while (++fit != fend);


    // adjust end iterators
    ++vend; ++fend;


    // collapse new edges
    vit=vend; vend=P.vertices_end();
    for (; vit!=vend; ++vit)
      halfedge_collapse(vit->halfedge()->opposite());
  }


  void halfedge_collapse(Polyhedron::Halfedge_descriptor pq)
  {
    // this code is copied from the CGAL surface simplification package

    Polyhedron::Halfedge_descriptor qp = pq->opposite();
    Polyhedron::Halfedge_descriptor pt = pq->prev()->opposite();
    Polyhedron::Halfedge_descriptor qb = qp->prev()->opposite();

    bool lTopFaceExists         = !pq->is_border() ;
    bool lBottomFaceExists      = !qp->is_border() ;
    bool lTopLeftFaceExists     = lTopFaceExists    && !pt->is_border() ;
    bool lBottomRightFaceExists = lBottomFaceExists && !qb->is_border() ;

    Polyhedron::Vertex_descriptor q = pq->vertex();
    Polyhedron::Vertex_descriptor p = pq->opposite()->vertex();

    bool lP_Erased = false, lQ_Erased = false ;

    if ( lTopFaceExists )
    {
      if ( lTopLeftFaceExists )
      {
        P.join_facet (pt);
      }
      else
      {
        P.erase_facet(pt->opposite());

        if ( !lBottomFaceExists )
        {
          lP_Erased = true ;
        }
      }
    }

    if ( lBottomFaceExists )
    {
      if ( lBottomRightFaceExists )
      {
        P.join_facet (qb);
      }
      else
      {
        P.erase_facet(qb->opposite());

        if ( !lTopFaceExists )
        {
          lQ_Erased = true ;
        }
      }
    }

    if ( !lP_Erased && !lQ_Erased )
    {
      P.join_vertex(pq);
      lP_Erased = true ;
    }
  }

};


//=============================================================================
