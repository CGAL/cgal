/***************************************************************************
cgal_types.h  -  description
                             -------------------
begin                : jan 2002
copyright            : (C) 2002 by Pierre Alliez - INRIA
email                : pierre.alliez@sophia.inria.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <CGAL/basic.h>

// LS 12/2004: seems mandatory to include gl.h
#ifdef WIN32
#include <windows.h>
#endif

// CGAL stuff

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Color.h>
#include <algorithm>
#include <vector>
#include <list>
#include <GL/gl.h>
#include <stdio.h>

#include "my_kernel.h"
#include "Feature_skeleton.h"

#ifndef _CGAL_DEFS_
#define _CGAL_DEFS_

typedef My_kernel::Vector_3 Vector;
typedef My_kernel::Vector_2 Vector_2;
typedef My_kernel::Point_3  Point;
typedef My_kernel::Point_2  Point_2;
typedef My_kernel::Segment_2 Segment_2;

class My_sample
{
public:
  enum type_t {TYPE_CORNER,
               TYPE_FEATURE,
               TYPE_SURFACE};
  enum location_t {LOCATION_VERTEX,
                   LOCATION_BORDER,
                   LOCATION_INNER};

private:

  type_t m_type; // feature/surface
  location_t m_location; // inner/border

  // coordinates
  double m_coord[3];

  // in parametric space
  double m_u;
  double m_v;

  // for feature samples
  double m_alpha; // barycentric coordinates
  int m_index_halfedge;

  // twin sample
  My_sample *m_pTwin;

  // density
  double m_density_uv;

  // located facet
  void *m_pFacetLocate;

public:
  My_sample(const double x,
            const double y,
            const double z,
            type_t t,
            location_t l,
            double density_uv)
  {
    m_coord[0] = x;
    m_coord[1] = y;
    m_coord[2] = z;
    m_type = t;
    m_location = l;
    m_u = m_v = 0;
    m_index_halfedge = -1;
    m_alpha = 0.0;
    m_density_uv = density_uv;
    m_pFacetLocate = NULL;
  }

  My_sample(const My_sample &sample)
  {
    m_coord[0] = sample.x();
    m_coord[1] = sample.y();
    m_coord[2] = sample.z();
    m_type = sample.type();
    m_location = sample.location();
    m_u = sample.u();
    m_v = sample.v();
    m_index_halfedge = sample.index_halfedge();
    m_alpha = sample.alpha();
    m_density_uv = sample.density_uv();
    m_pFacetLocate = NULL;

    /*
    std::cerr << "type: " << sample.type() << std::endl;
    std::cerr << "location: " << sample.location() << std::endl;
    std::cerr << "uv: " << sample.u() << " " << sample.v() << std::endl;*/
  }

  My_sample(const double x,
            const double y,
            const double z,
            const double u,
            const double v,
            type_t type,
            location_t location,
            int index_halfedge,
            double alpha,
            double density_uv)
  {
    CGAL_assertion(type == TYPE_FEATURE);
    m_coord[0] = x;
    m_coord[1] = y;
    m_coord[2] = z;
    m_u = u;
    m_v = v;
    m_type = type;
    m_alpha = alpha;
    m_location = location;
    m_index_halfedge = index_halfedge;
    m_density_uv = density_uv;
    CGAL_assertion(m_index_halfedge >= 0);
    m_pFacetLocate = NULL;
  }

  My_sample(const double x,
            const double y,
            const double z,
            const double u,
            const double v,
            type_t type,
            location_t location,
            double density_uv)
  {
    m_coord[0] = x;
    m_coord[1] = y;
    m_coord[2] = z;
    m_u = u;
    m_v = v;
    m_type = type;
    m_alpha = 0.0;
    m_location = location;
    m_index_halfedge = -1;
    m_density_uv = density_uv;
    m_pFacetLocate = NULL;
  }

  // twin
  My_sample *twin() { return m_pTwin; }
  //const My_sample *twin() const { return m_pTwin; }
  void twin(My_sample *pTwin) { m_pTwin = pTwin; }

  // located facet
  void facet_locate(void *pFacet) { m_pFacetLocate = pFacet; }
  void *facet_locate() { return m_pFacetLocate; }

  // coordinates
  double x() { return m_coord[0]; }
  double y() { return m_coord[1]; }
  double z() { return m_coord[2]; }
  const double x() const { return m_coord[0]; }
  const double y() const { return m_coord[1]; }
  const double z() const { return m_coord[2]; }
  void x(double x) { m_coord[0] = x; }
  void y(double y) { m_coord[1] = y; }
  void z(double z) { m_coord[2] = z; }
  void xyz(const double x,
           const double y,
           const double z)

  {
    m_coord[0] = x;
    m_coord[1] = y;
    m_coord[2] = z;
  }
  void xyz(const double coord[])

  {
    m_coord[0] = coord[0];
    m_coord[1] = coord[1];
    m_coord[2] = coord[2];
  }
  double operator[](unsigned int i)
  {
    CGAL_assertion(i <= 2);
    return m_coord[i];
  }

  // parametric coordinates
  double u() { return m_u; }
  double v() { return m_v; }
  const double u() const { return m_u; }
  const double v() const { return m_v; }
  void u(double u) { m_u = u; }
  void v(double v) { m_v = v; }
  void uv(double u,double v)
  {
    m_u = u;
    m_v = v;
  }

  // density
  double density_uv() { return m_density_uv; }
  const double density_uv() const { return m_density_uv; }
  void density_uv(double density_uv) { m_density_uv = density_uv; }

  // for univariate lloyd
  double alpha() { return m_alpha; }
  const double alpha() const { return m_alpha; }
  void alpha(const double alpha) { m_alpha = alpha; }

  int index_halfedge() { return m_index_halfedge; }
  const int index_halfedge() const { return m_index_halfedge; }
  void index_halfedge(const int index_halfedge) { m_index_halfedge = index_halfedge; }

  // location and type
  location_t location() { return m_location; }
  type_t type() { return m_type; }
  const location_t location() const { return m_location; }
  const type_t type() const { return m_type; }
};

// Two functors to compute the normals:  We assume the
// Simple_cartesian<double> My_kernel here and use its global functions.

// compute facet normal
struct Facet_normal
{
  template <class Facet>
  void operator()(Facet& f)
  {
    typename Facet::Normal_3 sum = CGAL::NULL_VECTOR;
    typename Facet::Halfedge_around_facet_circulator h = f.facet_begin();
    do
    {
      typename Facet::Normal_3 normal = CGAL::cross_product(
        h->next()->vertex()->point() - h->vertex()->point(),
        h->next()->next()->vertex()->point() - h->next()->vertex()->point());
      double sqnorm = normal * normal;
      if(sqnorm != 0)
        normal = normal / std::sqrt(sqnorm);
      sum = sum + normal;
    }
    while(++h != f.facet_begin());
    double sqnorm = sum * sum;
    if(sqnorm != 0.0)
      f.normal() = sum / std::sqrt(sqnorm);
    else
      std::cerr << "degenerate face" << std::endl;
    // NULL << f.normal().x() << " " <<
    //         f.normal().y() << " " <<
    //         f.normal().z() << std::endl;
  }
};

// compute facet center
struct Facet_center
{
    template <class Facet>
    void operator()(Facet& f)
    {
        Vector vec(0.0,0.0,0.0);
        int degree = 0;
        typedef typename Facet::Halfedge_around_facet_const_circulator circ;
        circ h = f.facet_begin();
        do
        {
          vec = vec + (h->vertex()->point()-CGAL::ORIGIN);
          degree++;
        }
        while (++h != f.facet_begin());
        f.center() = CGAL::ORIGIN + (vec/degree);
    }
};

// compute facet area
struct Facet_area
{
    template <class Facet>
    void operator()(Facet& f)
    {
        typedef typename Facet::Halfedge_around_facet_const_circulator circ;
        circ h = f.facet_begin();
        
        // check min. degree
        int degree = 0;
        do
          degree++;
        while(++h != f.facet_begin());
        CGAL_assertion(degree >= 3);

        // ** TODO ** for polygons
        Point p0 = h->vertex()->point();
        Point p1 = h->next()->vertex()->point();
        Point p2 = h->next()->next()->vertex()->point();
        Vector bc = p2-p1;
        Vector ba = p0-p1;
        Vector v = CGAL::cross_product(bc,ba);
        f.area(0.5 * std::sqrt(v*v));
    }
};


struct Vertex_normal
{
    template <class Vertex>
    void operator()( Vertex& v)
    {
        typename Vertex::Normal_3 normal = CGAL::NULL_VECTOR;
        typedef typename Vertex::Halfedge_around_vertex_const_circulator Circ;
        Circ c = v.vertex_begin();
        Circ d = c;
        CGAL_For_all( c, d) {
            if ( ! c->is_border())
                normal = normal + c->facet()->normal();
        }
        v.normal() = normal / std::sqrt( normal * normal);
    }
};


template <class Refs, class T, class Norm>
class My_facet : public CGAL::HalfedgeDS_face_base<Refs, T>
{
  // face data
  Norm m_normal;
  CGAL::Color m_color;
  int m_tag;
  Point m_center;

  // information for sampling
  double m_area;
  double m_geometry; // total amount of geometry
  double m_geometry_corrected;

  // param
  double m_area_uv;    // area in parametric space
  double m_stretching; // area stretching

public:

  // life cycle
  // no constructors to repeat, since only
  // default constructor mandatory

  My_facet()
  {
    m_area = 0.0;
    m_geometry = 0.0;
    m_color = CGAL::WHITE;
    m_tag = 0;
  }

  // normal
  typedef Norm Normal_3;
  Normal_3&       normal()       { return m_normal; }
  const Normal_3& normal() const { return m_normal; }

  // center
  Point& center() { return m_center; }
  const Point& center() const { return m_center; }

  // area
  double area() { return m_area; }
  void area(double a) { m_area = a; }
  const double area() const { return m_area; }

  // area_uv
  double area_uv() { return m_area_uv; }
  void area_uv(double a) { m_area_uv = a; }
  const double area_uv() const { return m_area_uv; }

  // stretching
  double stretching() { return m_stretching; }
  void stretching(double s) { m_stretching = s; }

  // geometry
  double geometry() { return m_geometry; }
  void geometry(double g) { m_geometry = g; }
  void add_geometry(double x) { m_geometry_corrected += x; }
  const double geometry() const { return m_geometry; }

  double geometry_corrected() { return m_geometry_corrected; }
  const double geometry_corrected() const { return m_geometry_corrected; }
  void geometry_corrected(double g) { m_geometry_corrected = g; }

  // color
  void color(CGAL::Color &color) { m_color = color; }
  CGAL::Color& color() { return m_color; }
  const CGAL::Color& color() const { return m_color; }
  void color(unsigned char r,unsigned char g,unsigned char b)
  {
    m_color = CGAL::Color(r,g,b);
  }

  // tag
  int tag() { return m_tag; }
  const int tag() const { return m_tag; }
  void tag(int tag) { m_tag = tag; }

  // distance
  double distance(Point *pPoint)
  {
    Vector vec = (*pPoint-m_center);
    return My_kernel::len(vec);
  }
};

template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class My_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
  Norm m_normal;
  bool m_is_sharp;
  int m_tag;

  // information for sampling
  double m_len;
  double m_geometry; // total amount of geometry
  double m_geometry_corrected;

  // parameterization
  bool m_is_parameterized;
  bool m_is_seaming; // seaming halfedge
  // texture coordinates
  double m_u;
  double m_v;
  int m_index; // for param.

  // surface cutting
  float m_distance;
  char m_superimpose;

public:

  // life cycle
  // no constructors to repeat, since only
  // default constructor mandatory

  My_halfedge()
  {
    m_tag = 0;
    m_u = 0.0;
    m_v = 0.0;
    m_len = 0.0;
    m_index = 0;
    m_geometry = 0.0;
    m_is_sharp = false;
    m_is_seaming = false;
    m_is_parameterized = false;
    m_superimpose = false;
  }

  // tag
  int tag() { return m_tag; }
  const int tag() const { return m_tag; }
  void tag(int tag) { m_tag = tag; }

  // sharpness
  bool is_sharp() { return m_is_sharp; }
  void is_sharp(bool sharp) { m_is_sharp = sharp; }

  // seaming backbone
  bool is_seaming() { return m_is_seaming; }
  void is_seaming(bool is_seaming) { m_is_seaming = is_seaming; }

  // edge superimposing
  char superimpose() { return m_superimpose; }
  void superimpose(char s) { m_superimpose = s; }

  // precomputed distance
  float distance() { return m_distance; }
  void distance(float distance) { m_distance = distance; }

  // normal
  typedef Norm Normal_3;
  Normal_3&       normal()       { return m_normal; }
  const Normal_3& normal() const { return m_normal; }

  // texture coordinates
  double u() { return m_u; }
  double v() { return m_v; }
  void u(double u) { m_u = u; }
  void v(double v) { m_v = v; }
  void uv(double u,double v) { m_u = u; m_v = v; }
  double len_uv()
  {
    double pu = prev()->u();
    double pv = prev()->v();
    return sqrt((pu-m_u)*(pu-m_u)+(pv-m_v)*(pv-m_v));
  }
  double len() { return m_len; }
  void len(double l) { m_len = l; }

  // param.
  bool is_parameterized()  { return m_is_parameterized; }
  void is_parameterized(bool is)  { m_is_parameterized = is; }

  // index
  int index() { return m_index; }
  void index(int i) { m_index = i; }

  // geometry
  double geometry() { return m_geometry; }
  void geometry(double g) { m_geometry = g; }
  void add_geometry(double x) { m_geometry_corrected += x; }
  const double geometry() const { return m_geometry; }

  double geometry_corrected() { return m_geometry_corrected; }
  void geometry_corrected(double g) { m_geometry_corrected = g; }
  const double geometry_corrected() const { return m_geometry_corrected; }
};


// A redefined items class for the Polyhedron_3 with a refined vertex
// class that contains a member for the normal vector and a refined
// facet with a normal vector instead of the plane equation (this is
// an alternative solution instead of using Polyhedron_traits_with_normals_3).

template <class Refs, class T, class P, class Norm>
class My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
  // normal
  Norm m_normal;

  // texture coordinates
  double m_u;
  double m_v;

  // multiplicity on the seaming backbone
  unsigned char m_multiplicity;
  // the vertex can be on the boundary or inner
  bool m_is_inner_param;

  // mean curvature normal & density
  Vector m_meanCurvatureNormal;
  double m_meanCurvature;
  double m_density;
  double m_density_uv;
  double m_area_around;

  // feature skeleton
  int m_sharp_degree;
  bool m_is_crease;
  bool m_is_corner;
  bool m_is_border;

  // stretching factor on incident facets
  double m_stretching;

  // index
  int m_index;

  // misc
  CGAL::Color m_color;
  int m_tag;

  // anisotropic data (one tensor per vertex)
  double m_kmin;
  double m_kmax;
  double m_orientation_kmax;

  // reference to direction system
  void *m_pDirection;

public:
  // life cycle
  My_vertex()  { init(); }
  // repeat mandatory constructors
  My_vertex(const P& pt)
    : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt)
  {
    init();
  }

  // direction
  void *direction() { return m_pDirection; }
  void direction(void *pDirection) { m_pDirection = pDirection; }

  // anisotropic data
  double kmin() { return m_kmin; }
  double kmax() { return m_kmax; }
  double orientation_kmax() { return m_orientation_kmax; }
  void kmin(double kmin) { m_kmin = kmin; }
  void kmax(double kmax) { m_kmax = kmax; }
  void orientation_kmax(double o) { m_orientation_kmax = o; }
  void area_around(double area) { m_area_around = area; }
  double area_around() { return m_area_around; }

  // stretching
  void stretching(double stretching) { m_stretching = stretching; }
  double stretching() { return m_stretching; }

  void init()
  {
    m_meanCurvatureNormal = CGAL::NULL_VECTOR;
    m_meanCurvature = 1.0; // unitary
    m_is_corner = false;
    m_is_border = false;
    m_is_crease = false;
    m_sharp_degree = 0;
    m_u = 0.0;
    m_v = 0.0;
    m_multiplicity = 1;
    m_is_inner_param = true;
    m_multiplicity = 0;
    m_pDirection = 0;
  }

  // sharp degree / corner
  void is_corner(bool is) { m_is_corner = is; }
  bool is_corner() { return m_is_corner; }
  void is_crease(bool is) { m_is_crease = is; }
  bool is_crease() { return m_is_crease; }
  bool is_feature() { return m_is_crease || m_is_corner; }
  void sharp_degree(int d) { m_sharp_degree = d; }
  int sharp_degree() { return m_sharp_degree; }

  // inner vertex for the parametrization
  bool is_inner_param() { return m_is_inner_param; }
  void is_inner_param(bool inner) { m_is_inner_param = inner; }

  // multiplicity
  unsigned char multiplicity() { return m_multiplicity; }
  void multiplicity(unsigned char m) { m_multiplicity = m; }
  void increase_multiplicity() { m_multiplicity++; }

  // index
  int index() { return m_index; }
  void index(int i) { m_index = i; }

  // texture coordinates
  double u() { return m_u; }
  double v() { return m_v; }
  void u(double u) { m_u = u; }
  void v(double v) { m_v = v; }
  void uv(double u,double v) { m_u = u; m_v = v; }

  // normal
  typedef Norm Normal_3;
  Normal_3&       normal()       { return m_normal; }
  const Normal_3& normal() const { return m_normal; }

  // mean curvature normal
  Vector&       meanCurvatureNormal()       { return m_meanCurvatureNormal; }
  const Vector& meanCurvatureNormal() const { return m_meanCurvatureNormal; }

  // density
  double density() { return m_density; }
  void density(double density) { m_density = density; }
  double density_uv() { return m_density_uv; }
  void density_uv(double density_uv) { m_density_uv = density_uv; }

  // color
  CGAL::Color& color() { return m_color; }
  const CGAL::Color& color() const { return m_color; }

  // tag
  int tag() { return m_tag; }
  const int tag() const { return m_tag; }
  void tag(int tag) { m_tag = tag; }

  // curvature
  double meanCurvature() { return m_meanCurvature; }
  void meanCurvature(double mc) { m_meanCurvature = mc; }
  const double meanCurvature() const { return m_meanCurvature; }
};

struct My_items : public CGAL::Polyhedron_items_3
{
    // wrap vertex
    template <class Refs, class Traits>
    struct Vertex_wrapper
    {
        typedef typename Traits::Point_3  Point;
        typedef typename Traits::Vector_3 Normal;
        typedef My_vertex<Refs,
                          CGAL::Tag_true,
                          Point,
                          Normal> Vertex;
    };

    // wrap face
    template <class Refs, class Traits>
    struct Face_wrapper
    {
        typedef typename Traits::Vector_3 Normal;
        typedef My_facet<Refs,
                         CGAL::Tag_true,
                         Normal> Face;
    };

    // wrap halfedge
    template <class Refs, class Traits>
    struct Halfedge_wrapper
    {
        typedef typename Traits::Vector_3 Normal;
        typedef My_halfedge<Refs,
                            CGAL::Tag_true,
                            CGAL::Tag_true,
                            CGAL::Tag_true,
                            Normal> Halfedge;
    };
};


typedef CGAL::Polyhedron_3<My_kernel,My_items> Polyhedron;

class Polyhedron_ex : public Polyhedron
{
  private:

    // rendering
    //********************************
    int m_gl_list;
    int m_view;

    // remeshing stuffs
    double m_scaleEnvelope;
    std::list<My_sample> m_samples;

  public:

  /*
    typedef Polyhedron::Vertex                                   Vertex;
    typedef Polyhedron::Vertex_handle                            Vertex_handle;
    typedef Polyhedron::Vertex_iterator                          Vertex_iterator;
    typedef Polyhedron::Halfedge_handle                          Halfedge_handle;
    typedef Polyhedron::Edge_iterator                            Edge_iterator;
    typedef Polyhedron::Facet_iterator                           Facet_iterator;
    typedef Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;
    typedef Polyhedron::Halfedge_around_facet_circulator         HF_circulator;  */

 private:

    // feature/boundary skeletons
    typedef Feature_backbone<Vertex_handle,Halfedge_handle> backbone;
    typedef Feature_skeleton<Vertex_handle,Halfedge_handle> skeleton;
    skeleton m_skeleton;
    backbone m_seaming_backbone;
  
  public:

    // data access
    skeleton* get_skeleton() { return &m_skeleton; }
    backbone* get_seaming_backbone() { return &m_seaming_backbone; }
  
    // life cycle
    //********************************
    Polyhedron_ex()
    {
    }
    virtual ~Polyhedron_ex()
    {
      free_skeleton();
    }
    
    void free_skeleton()
    {
      m_samples.clear();
      m_skeleton.free();
      //m_seaming_backbone.halfedges()->clear();
    }
    
    // data access
    //********************************
    std::list<My_sample> *samples() { return &m_samples; }
    
    // normals
    //********************************
    void compute_normals()
    {
      fprintf(stderr,"  compute normals...");
      compute_normals_per_facet();
      compute_normals_per_vertex();
      fprintf(stderr,"ok\n");
    }
    void compute_normals_per_facet()
    {
      std::for_each(facets_begin(),facets_end(),Facet_normal());
    }
    void compute_normals_per_vertex()
    {
      std::for_each(vertices_begin(),vertices_end(),Vertex_normal());
    }
    void compute_facet_areas()
    {
      fprintf(stderr,"  compute facet areas...");
      std::for_each(facets_begin(),facets_end(),Facet_area());
      fprintf(stderr,"ok\n");
    }
    
    void compute_facet_uv_areas()
    {
      fprintf(stderr,"  compute uv facet areas...");
      Facet_iterator pFace;
      for(pFace = facets_begin();pFace != facets_end();pFace++)
      {
        Halfedge_around_facet_circulator h = pFace->facet_begin();

        // a -> [10]
        // b -> [20]
        // a.u*b.v-b.u*a.v
        double au = h->u()-h->next()->u();          
        double av = h->v()-h->next()->v();          
        double bu = h->u()-h->next()->next()->u();          
        double bv = h->v()-h->next()->next()->v();          
        double area = 0.5*(au*bv-bu*av);
        pFace->area_uv(area);
      }
      fprintf(stderr,"ok\n");
    }
    
    void compute_facet_stretching()
    {
      compute_facet_areas();
      compute_facet_uv_areas();
      Facet_iterator pFace;
      for(pFace = facets_begin();pFace != facets_end();pFace++)
      {
        double area = pFace->area();
        double area_uv = pFace->area_uv();
        if(area_uv != 0.0)
          pFace->stretching(area/area_uv);
        else
          pFace->stretching(0.0);
      }
    }
    
    void compute_vertex_stretching()
    {
      Vertex_iterator pVertex;
      for(pVertex = vertices_begin();
          pVertex != vertices_end();
          pVertex++)
      {
        double sum_2d = 0.0;
        double sum_3d = 0.0;
        typedef Halfedge_around_vertex_circulator circ;
        circ pHalfedge = pVertex->vertex_begin();
        circ end = pHalfedge;
        CGAL_For_all(pHalfedge,end)
        {
          if(!pHalfedge->is_border())
          {
            sum_2d += fabs(pHalfedge->facet()->area_uv());
            sum_3d += pHalfedge->facet()->area();
          }
        }
        if(sum_2d > 0.0)
          pVertex->stretching(sum_3d/sum_2d);
        else
          pVertex->stretching(1.0);
      }
    }
    
    // facet centers
    void compute_facet_centers()
    {
      fprintf(stderr,"  compute facet centers...");
      std::for_each(facets_begin(),facets_end(),Facet_center());
      fprintf(stderr,"ok\n");
    }

    // color all facets
    void color_facets(CGAL::Color &color)
    {
      std::cerr << "  color all facets...";
      Facet_iterator pFace;
      for(pFace = facets_begin();pFace != facets_end();pFace++)
        pFace->color(color);
      std::cerr << "ok" << std::endl;
    }
    
    // color all facets
    void color_facets(const unsigned char r,
                      const unsigned char g,
                      const unsigned char b)
    {
      std::cerr << "  color all facets...";
      Facet_iterator pFace;
      for(pFace = facets_begin();pFace != facets_end();pFace++)
        pFace->color(r,g,b);
      std::cerr << "ok" << std::endl;
    }
    
    // tag all facets
    void tag_facets(const int tag)
    {
      Facet_iterator pFace;
      for(pFace = facets_begin();
          pFace != facets_end();
          pFace++)
        pFace->tag(tag);
    }
    
    // get any inner facet
    Facet_handle get_any_inner_facet()
    {
      Facet_iterator pFace;
      for(pFace = facets_begin();
          pFace != facets_end();
          pFace++)
        if(is_inner(pFace))
          return pFace;
     return NULL;
    }

    // get any facet with tag
    Facet_handle get_any_facet_tag(int tag)
    {
      Facet_iterator pFace;
      for(pFace = facets_begin();
          pFace != facets_end();
          pFace++)
        if(pFace->tag() == tag)
          return pFace;
     return NULL;
    }

    // get closest inner face
    Facet_handle get_closest_inner_facet(Point *pPoint)
    {
      Facet_iterator pFace = facets_begin();
      Facet_handle pClosest = pFace;
      double min = pFace->distance(pPoint);
      for(;pFace != facets_end();
          pFace++)
      {
        if(is_inner(pFace))
        {
          double distance = pFace->distance(pPoint);
          if(distance < min)
          {
            pClosest = pFace;
            min = distance;
          }
        }
      }
      return pClosest;
    }

    // tag all vertices
    bool is_inner(Facet_handle pFace)
    {
       typedef Halfedge_around_facet_const_circulator circ;
       circ h = pFace->facet_begin();
      do
      {
        if(h->opposite()->is_border())
          return false;
      }
      while(++h != pFace->facet_begin());
      return true;
    }
    
    // tag all vertices
    void tag_vertices(const int tag)
    {
      Vertex_iterator iter;
      for(iter = vertices_begin();
          iter != vertices_end();
          iter++)
      {
        Vertex_handle hVertex = iter;
        hVertex->tag(tag);
      }
    }
    
    // init vertex multiplicity
    void init_vertex_multiplicity(unsigned char multiplicity)
    {
      Vertex_iterator pVertex;
      for(pVertex = vertices_begin();
          pVertex != vertices_end();
          pVertex++)
        pVertex->multiplicity(multiplicity);
    }
    
    // trace vertex multiplicity
    void trace_vertex_multiplicity()
    {
      unsigned char max = 0;
      Vertex_iterator pVertex;
      for(pVertex = vertices_begin();
          pVertex != vertices_end();
          pVertex++)
        max = (pVertex->multiplicity() > max) ? pVertex->multiplicity() : max;
      if(max == 0)
        return;

      // init      
      unsigned int *pDistribution = new unsigned int[max+1];
      for(unsigned char i=0;i<=max;i++)
        pDistribution[i] = 0;
      
      // compute
      for(pVertex = vertices_begin();
          pVertex != vertices_end();
          pVertex++)
        pDistribution[pVertex->multiplicity()]++;
      
      // trace
      std::cerr << "    vertex multiplicity: ";
      for(unsigned char i=0;i<=max;i++)
        std::cerr << pDistribution[i] << " ";
      std::cerr << std::endl;

      // cleanup      
      delete [] pDistribution;
    }
    
    // compute bounding interval
    double min(int coord)
    {
      CGAL_assertion(size_of_vertices() > 0);
      Vertex_iterator pVertex = vertices_begin();
      double min = pVertex->point()[coord];
      for(;pVertex != vertices_end();pVertex++)
        min = std::min(min,pVertex->point()[coord]);
     return min;
    }
    double max(int coord)
    {
      CGAL_assertion(size_of_vertices() > 0);
      Vertex_iterator pVertex = vertices_begin();
      double max = pVertex->point()[coord];
      for(;pVertex != vertices_end();pVertex++)
        max = std::max(max,pVertex->point()[coord]);
     return max;
    }
    Vertex_handle vertex_min(int coord,
                             double &min)
    {
      CGAL_assertion(size_of_vertices() > 0);
      Vertex_iterator pVertex = vertices_begin();
      Vertex_handle pBest = pVertex;
      min = pVertex->point()[coord];
      for(;pVertex != vertices_end();pVertex++)
      {
        double value = pVertex->point()[coord];
        if(value < min)
        {
          min = std::min(min,value);
          pBest = pVertex;
        }
      }
     return pBest;
    }
    Vertex_handle vertex_max(int coord,
                             double &max)
    {
      CGAL_assertion(size_of_vertices() > 0);
      Vertex_iterator pVertex = vertices_begin();
      Vertex_handle pBest = pVertex;
      max = pVertex->point()[coord];
      for(;pVertex != vertices_end();pVertex++)
      {
        double value = pVertex->point()[coord];
        if(value > max)
        {
          max = std::max(max,value);
          pBest = pVertex;
        }
      }
     return pBest;
    }
    
    void set_is_inner_param_vertices(bool is_inner)
    {
      Vertex_iterator pVertex;
      for(pVertex = vertices_begin();
          pVertex != vertices_end();
          pVertex++)
        pVertex->is_inner_param(is_inner);
    }
    
    // tag all halfedges
    void tag_halfedges(const int tag)
    {
      Halfedge_iterator pHalfedge;
      for(pHalfedge = halfedges_begin();
          pHalfedge != halfedges_end();
          pHalfedge++)
        pHalfedge->tag(tag);
    }

    // superimposing tagging for all halfedges
    void superimpose_halfedges(const char s)
    {
      Halfedge_iterator pHalfedge;
      for(pHalfedge = halfedges_begin();
          pHalfedge != halfedges_end();
          pHalfedge++)
        pHalfedge->superimpose(s);
    }
    
    // flag all halfedges
    void flag_halfedges_seaming(const bool flag)
    {
      Halfedge_iterator pHalfedge;
      for(pHalfedge = halfedges_begin();
          pHalfedge != halfedges_end();
          pHalfedge++)
        pHalfedge->is_seaming(flag);
    }

    // valence of a vertex
    static unsigned int valence(Vertex_handle pVertex)
    {
      Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
      Halfedge_around_vertex_circulator end = pHalfEdge;
      unsigned int valence = 0;
      CGAL_For_all(pHalfEdge,end)
        valence++;
      return valence;
    }

    // degree of a face
    static unsigned int degree(Facet_handle pFace)
    {
      Halfedge_around_facet_circulator pHalfEdge = pFace->facet_begin();
      Halfedge_around_facet_circulator end = pHalfEdge;
      unsigned int valence = 0;
      CGAL_For_all(pHalfEdge,end)
        valence++;
      return valence;
    }

    // set all halfedges to is_param
    void set_is_parameterized_halfedges(const bool is)
    {
      Halfedge_iterator pHalfedge;
      for(pHalfedge = halfedges_begin();
          pHalfedge != halfedges_end();
          pHalfedge++)
        pHalfedge->is_parameterized(is);
    }
    
    int nb_adj_facets_tag(Vertex_handle pVertex,
                          int tag)
    {
      int nb = 0;
      typedef Halfedge_around_vertex_circulator circ;
      circ c = pVertex->vertex_begin();
      circ d = c;
      CGAL_For_all(c,d)
      {
        Halfedge_handle he = c;
        Facet_handle pFacet = he->facet();
        if(pFacet->tag() == tag)
          nb++;
      }
      return nb;
    }
    
    int nb_inner_param_vertices()
    {
      int nb = 0;
      Vertex_iterator pVertex;
      for(pVertex = vertices_begin();
          pVertex != vertices_end();
          pVertex++)
        if(pVertex->is_inner_param() == true)
          nb++;
      return nb;
    }
    
    // Number all mesh vertices following the order of the vertices_begin() iterator
    void precompute_vertex_indices()
    {
      Vertex_iterator pVertex;
      unsigned int i = 0;
      for(pVertex = vertices_begin();
          pVertex != vertices_end();
          pVertex++)
        pVertex->index(i++);
    }
    
    // Number all mesh half edges following the order of the halfedges_begin() iterator
    void precompute_halfedge_indices()
    {
      Halfedge_iterator pHalfedge;
      unsigned int i = 0;
      for(pHalfedge = halfedges_begin();
          pHalfedge != halfedges_end();
          pHalfedge++)
        pHalfedge->index(i++);
    }
    
    
    ///////////////////////////////////////////////////////////
    // GEOMETRY
    ///////////////////////////////////////////////////////////
    
    void operator*=(double scale)
    {
      std::cerr << "  scale..." << scale;
      Vertex_iterator pVertex;
      for(pVertex = vertices_begin();
          pVertex != vertices_end();
          pVertex++)
      {
        Vector v = pVertex->point() - CGAL::ORIGIN;
        v = v * scale;
        pVertex->point() = CGAL::ORIGIN + v;
      }
      std::cerr << "...ok" << std::endl;
    }
    void shift(double dx,
               double dy,
               double dz)
    {
      Vertex_iterator pVertex;
      for(pVertex = vertices_begin();
          pVertex != vertices_end();
          pVertex++)
        pVertex->point() = Point(pVertex->point().x()+dx,
                                 pVertex->point().y()+dy,
                                 pVertex->point().z()+dz);
    }
    void scale(double scale) { operator*=(scale); }

    // compute sharp edges wrt the dihedral angle
    void compute_sharp_edges(double dihedral_angle)
    {
      std::cerr << "  compute sharp edges...";
      Halfedge_iterator pHalfedge;
      int counter = 0;
      for(pHalfedge = halfedges_begin();
          pHalfedge != halfedges_end();
          pHalfedge++)
      {
        // border halfedges
        if(pHalfedge->facet() == NULL ||
           pHalfedge->opposite()->facet() == NULL)
          continue;

        Vector n1 = pHalfedge->facet()->normal();
        Vector n2 = pHalfedge->opposite()->facet()->normal();
        double angle = My_kernel::evaluate_angle(n1,n2);
        bool sharp = (angle > dihedral_angle);
        pHalfedge->is_sharp(sharp);
        if(sharp)
          counter++;
      }
      std::cerr << "..." << counter << "/" << size_of_halfedges() << "...ok" << std::endl;
    }

    // compute sharp edges wrt the dihedral angle
    void reset_sharp_edges(bool set)
    {
      std::cerr << "  reset sharp edges.";
      Halfedge_iterator pHalfedge;
      for(pHalfedge = halfedges_begin();
          pHalfedge != halfedges_end();
          pHalfedge++)
        pHalfedge->is_sharp(set);
      std::cerr << "..ok" << std::endl;
    }

    // compute corners (requires flags set on halfedges)
    unsigned int compute_corners(bool include_multiple_vertices,
                                 double angle)
    {
      std::cerr << "  compute...";
      unsigned int nb_corners = 0;
      unsigned int nb_creases = 0;
      // compute displacements
      Polyhedron_ex::Vertex_iterator pVertex;
      for(pVertex = vertices_begin();
         (pVertex != vertices_end());
          pVertex++)
      {
        Polyhedron_ex::Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
        Polyhedron_ex::Halfedge_around_vertex_circulator d = pHalfEdge;
        unsigned int sharp_degree = 0;
        CGAL_For_all(pHalfEdge,d)
          if(pHalfEdge->is_sharp())
            sharp_degree ++;
        pVertex->sharp_degree(sharp_degree);
        bool is_corner = (sharp_degree == 1 || sharp_degree >= 3);
        pVertex->is_corner(is_corner);
        pVertex->is_crease(sharp_degree == 2);
        nb_corners += is_corner;
        nb_creases += (sharp_degree == 2);

        // check if multiple
        if(!is_corner && include_multiple_vertices)
        {
          unsigned char m = pVertex->multiplicity();
          if(m == 1 || m > 2)
          {
            is_corner = true;
            pVertex->is_corner(true);
            pVertex->is_crease(false);
            nb_corners++;
          }

          // big ambiguity
          // one vertex is a real corner if at least
          // one sharp edge does not coincide with a seaming halfedge edge
          if(m == 2 && sharp_degree >= 1)
          {
            Polyhedron_ex::Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
            Polyhedron_ex::Halfedge_around_vertex_circulator d = pHalfEdge;
            CGAL_For_all(pHalfEdge,d)
              if(pHalfEdge->is_sharp() &&
                 !pHalfEdge->is_seaming())
              {
                is_corner = true;
                pVertex->is_corner(true);
                pVertex->is_crease(false);
                nb_corners++;
              }
          }
        }

        // check if the crease case is curvy enough
        // to become a corner
        if(!is_corner && sharp_degree==2)
        {
          Polyhedron_ex::Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
          Polyhedron_ex::Halfedge_around_vertex_circulator d = pHalfEdge;
          Polyhedron_ex::Halfedge_handle he[2];
          int i=0;
          CGAL_For_all(pHalfEdge,d)
            if(pHalfEdge->is_sharp())
              he[i++] = pHalfEdge;
          // compute dot product
          Vector u1 = he[0]->prev()->vertex()->point()-pVertex->point();
          Vector u2 = pVertex->point()-he[1]->prev()->vertex()->point();
          double dot = u1*u2;
          if(dot<0)
          {
            pVertex->is_corner(true);
            pVertex->is_crease(false);
            nb_corners++;
          }
        }

        // one boundary vertex can be a corner
        if(!is_corner && is_border(*pVertex))
        {
          Polyhedron_ex::Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
          Polyhedron_ex::Halfedge_around_vertex_circulator d = pHalfEdge;
          Polyhedron_ex::Halfedge_handle he[2];
          int i=0;
          CGAL_For_all(pHalfEdge,d)
            if(pHalfEdge->is_border() ||
               pHalfEdge->opposite()->is_border())
              he[i++] = pHalfEdge;

          // compute angle
          Vector u1 = he[0]->prev()->vertex()->point()-pVertex->point();
          Vector u2 = pVertex->point()-he[1]->prev()->vertex()->point();

          if(My_kernel::evaluate_angle(u1,u2) > angle)
          {
            pVertex->is_corner(true);
            pVertex->is_crease(false);
            nb_corners++;
          }
        }
      }
      std::cerr << "..." << nb_corners << " corners";
      std::cerr << "..." << nb_creases << " creases...ok" << std::endl;
      return nb_corners;
    }

    // propagate mean curvature from inner parts to
    // boundaries
    void propagate_mean_curvature()
    {
      Vertex_iterator pVertex;
      for(pVertex  = vertices_begin();
          pVertex != vertices_end();
          pVertex++)
      {
        if(is_border(*pVertex))
        {
          Halfedge_around_vertex_circulator hec = pVertex->vertex_begin();
          Halfedge_around_vertex_circulator d = hec;
          double max_curvature = 0.0;
          CGAL_For_all(hec, d)
          {
            Vertex_handle pVertex = hec->opposite()->vertex();
            if(!is_border(*pVertex))
            {
              double curvature = pVertex->meanCurvature();
              if(curvature > max_curvature)
                max_curvature = curvature;
            }
          }
          if(pVertex->meanCurvature() < max_curvature)
            pVertex->meanCurvature(max_curvature);
        }
      }
    }
    
    void compute_mean_curvature_normal(bool normalize = true)
    {
      fprintf(stderr,"  compute mean curvature...");
      Vertex_iterator last_v = vertices_end();
      last_v--;
      Vertex_iterator v = vertices_begin();
      double max_curvature = 0.0;
      int counter_pure_voronoi = 0;
      do
      {
        bool pure_voronoi = 0;
        double curvature = compute_mean_curvature_normal(*v,pure_voronoi);
        counter_pure_voronoi += (pure_voronoi == true);
        max_curvature = std::max(max_curvature,curvature);
      }
      while(v++ != last_v);
      //fprintf(stderr,"ok (max: %g,pure voronoi: %d out of %d)\n",
      fprintf(stderr,"ok\n");
    
      // propagate mean curvature from inner parts to boundaries
      propagate_mean_curvature();
    
      // normalize
      if(normalize && max_curvature!= 0.0)
      {
        fprintf(stderr,"  normalize mean curvature...");
        v = vertices_begin();
        do
          v->meanCurvature(v->meanCurvature()/max_curvature);
        while(v++ != last_v);
        fprintf(stderr,"ok\n");
      }
    }

    // total area
    //********************************
    double total_area()
    {
      Facet_iterator last_f = facets_end();
      last_f--;
      Facet_iterator f = facets_begin();
      double total = 0.0;
      do
        total += f->area();
      while(f++ != last_f);
      return total;
    }

    // dump the param to an eps file
    //********************************
    bool dump_param(char *pFilename)
    {
      std::cerr << "dump parameterization to " << pFilename << "..." << std::endl;
      FILE *pFile = fopen(pFilename,"wt");
      if(pFile == NULL)
      {
        std::cerr << "  unable to open file " << pFilename <<  " for writing" << std::endl;
        return false;
      }
      dump_param(pFile);
      fclose(pFile);
      return true;
    }

    // dump the param to the stdout
    //********************************
    void dump_param()  { dump_param(stdout); }

    // dump the param to an eps file
    //********************************
    void dump_param(FILE *pFile,
                    double scale = 500.0)
    {
      // compute bounding box
      double xmin,xmax,ymin,ymax;
      xmin = ymin = xmax = ymax = 0;
      Halfedge_iterator pHalfedge;
      for(pHalfedge = halfedges_begin();
          pHalfedge != halfedges_end();
          pHalfedge++)
      {
        double x1 = scale * pHalfedge->prev()->u();
        double y1 = scale * pHalfedge->prev()->v();
        double x2 = scale * pHalfedge->u();
        double y2 = scale * pHalfedge->v();
        xmin = std::min(xmin,x1);
        xmin = std::min(xmin,x2);
        xmax = std::max(xmax,x1);
        xmax = std::max(xmax,x2);
        ymax = std::max(ymax,y1);
        ymax = std::max(ymax,y2);
        ymin = std::min(ymin,y1);
        ymin = std::min(ymin,y2);
      }

      fprintf(pFile,"%%!PS-Adobe-2.0 EPSF-2.0\n");
      fprintf(pFile,"%%%%BoundingBox: %g %g %g %g\n",xmin,ymin,xmax,ymax);
      fprintf(pFile,"%%%%EndComments\n");
      fprintf(pFile,"gsave\n");
      fprintf(pFile,"0.1 setlinewidth\n");

      // color macros
      fprintf(pFile,"\n%% RGB color command - r g b C\n");
      fprintf(pFile,"/C { setrgbcolor } bind def\n");
      fprintf(pFile,"/white { 1 1 1 C } bind def\n");
      fprintf(pFile,"/black { 0 0 0 C } bind def\n");

      // edge macro -> E
      fprintf(pFile,"\n%% Black stroke - x1 y1 x2 y2 E\n");
      fprintf(pFile,"/E {moveto lineto stroke} bind def\n");
      fprintf(pFile,"black\n\n");

      // for each halfedge
      for(pHalfedge = halfedges_begin();
          pHalfedge != halfedges_end();
          pHalfedge++)
      {
        double x1 = scale * pHalfedge->prev()->u();
        double y1 = scale * pHalfedge->prev()->v();
        double x2 = scale * pHalfedge->u();
        double y2 = scale * pHalfedge->v();
        fprintf(pFile,"%g %g %g %g E\n",x1,y1,x2,y2);
      }

      /* Emit EPS trailer. */
      fputs("grestore\n\n",pFile);
      fputs("showpage\n",pFile);
    }

    // output to an obj file
    //********************************
    bool write_file_obj(const char *pFilename)
    {
      std::cerr << "dump mesh to " << pFilename << "..." << std::endl;
      FILE *pFile = fopen(pFilename,"wt");
      if(pFile == NULL)
      {
        std::cerr << "  unable to open file " << pFilename <<  " for writing" << std::endl;
        return false;
      }
      
      bool ok = write_file_obj(pFile);

      fclose(pFile);
      
      return ok;
    }

	// output to an obj file
	// v x y z
	// f 1 2 3 4 (1-based)
	//********************************
	bool write_file_obj(FILE *pFile)
	{
		fprintf(stderr,"  write_file_obj()...");

		// Number all mesh vertices following the order of the vertices_begin() iterator
		precompute_vertex_indices();
		// Number all mesh half edges following the order of the halfedges_begin() iterator
		precompute_halfedge_indices();

		// write the name of material file
		fprintf(pFile, "mtllib parameterization.mtl\n") ;

		// output coordinates
		fprintf(pFile, "# vertices\n") ;
		Vertex_iterator pVertex;
		for(pVertex = vertices_begin();
			pVertex != vertices_end();
			pVertex++)
			fprintf(pFile,"v %g %g %g\n",
				pVertex->point().x(),
				pVertex->point().y(),
				pVertex->point().z());

		// Write UVs (1 UV / half edge)
		fprintf(pFile, "# uv coordinates\n") ;
		Halfedge_iterator pHalfedge;
		for(pHalfedge = halfedges_begin();
			pHalfedge != halfedges_end();
			pHalfedge++)
			fprintf(pFile,"vt %5.2f %5.2f\n",
					pHalfedge->u(),
					pHalfedge->v());

		// Write faces using the unique material # 1
		fprintf(pFile, "# faces\nusemtl Mat_1\n");
			Facet_iterator pFacet;
		for(pFacet = facets_begin();
			pFacet != facets_end();
			pFacet++)
		{
			Halfedge_around_facet_circulator h = pFacet->facet_begin();
			fprintf(pFile,"f ");
			do
				fprintf(pFile, "%d/%d ", h->vertex()->index()+1, h->index()+1);
			while(++h != pFacet->facet_begin());
			fprintf(pFile,"\n");
		}

		fprintf(stderr,"ok\n");

		return true;
    }

    // dump .obj file to stdout
    //********************************
    bool write_file_obj()  { return write_file_obj(stdout); }

  // boundary
  static bool is_border(const Vertex &v)
  {
    Halfedge_around_vertex_const_circulator c = v.vertex_begin();
    if(c == NULL) // isolated vertex
      return 1;
    Halfedge_around_vertex_const_circulator d = c;
    CGAL_For_all(c, d)
      if(c->is_border())
        return 1;
    return 0;
  }
  // is vertex on border ?
  static bool is_border(Vertex_handle pVertex)
  {
    Halfedge_around_vertex_const_circulator pHalfedge = pVertex->vertex_begin();
    Halfedge_around_vertex_const_circulator end = pHalfedge;
    if(pHalfedge == NULL) // isolated vertex
      return true;
    CGAL_For_all(pHalfedge,end)
      if(pHalfedge->is_border())
        return true;
    return false;
  }

  int find_index_first_parameterized_halfedge(Halfedge_handle seed,
                                              Vertex_handle pivot,
                                              Halfedge_handle &which)
  {
    Halfedge_handle halfedge = seed;
    while(1)
    {
      if(halfedge->is_parameterized())
      {
        which = halfedge;
        return halfedge->index();
      }
      // ** TODO ** triangle case only
      Halfedge_handle opposite = halfedge->opposite();
      if(opposite->is_border())
      {
        which = opposite->prev();
        CGAL_assertion(which->is_parameterized());
        return which->index();
      }
      else
        halfedge = halfedge->opposite()->next()->next();
    }
    return -1;
  }
  int find_index_first_halfedge_tag(Halfedge_handle seed,
                                    Vertex_handle pivot,
                                    Halfedge_handle &which,
                                    int tag)
  {
    Halfedge_handle halfedge = seed;
    while(1)
    {
      if(halfedge->tag() == tag)
      {
        which = halfedge;
        return halfedge->index();
      }
      Halfedge_handle opposite = halfedge->opposite();
      if(opposite->is_border())
      {
        which = opposite->prev();
        CGAL_assertion(which->tag() == tag);
        return which->index();
      }
      else
        halfedge = halfedge->opposite()->next()->next();
    }
    return -1;
  }

  // halfedge len
  double len(Polyhedron::Halfedge_handle halfedge)
  {
    Vector v = (halfedge->vertex()->point()-
                halfedge->prev()->vertex()->point());
    return std::sqrt(v*v);
  }

  // mean curvature
  // from http://www.multires.caltech.edu/pubs/diffGeoOps.pdf
  double compute_mean_curvature_normal(Vertex &v,
                                       bool &pure_voronoi)
  {
    // mean curvature normal on a boundary vertex
    Vector vec = CGAL::NULL_VECTOR;
    if(is_border(v))
    {
      // gather 2 adjacent boundary vertices
      Vertex_handle pVertices[2];
      Halfedge_around_vertex_circulator hec = v.vertex_begin();
      Halfedge_around_vertex_circulator d = hec;
      CGAL_assertion(hec != NULL);
      int nb = 0;
      CGAL_For_all(hec, d)
      {
        CGAL_assertion(hec != NULL);
        Vertex_handle pVertex = hec->opposite()->vertex();
        if(is_border(*pVertex) &&
           (hec->is_border() || hec->opposite()->is_border()))
          pVertices[nb++] = pVertex;
      }
      CGAL_assertion(nb == 2);

      // compute mean curvature normal (v is b)
      Vector ba = pVertices[0]->point()-v.point();
      Vector bc = pVertices[1]->point()-v.point();
      double norm_ba = My_kernel::len(ba);
      double norm_bc = My_kernel::len(bc);
      double sum = norm_ba+norm_bc;
      if(norm_ba != 0 && norm_bc != 0)
        vec = 0.5*(ba/norm_ba+bc/norm_bc)/sum;
      double curvature = My_kernel::len(vec);
      v.meanCurvatureNormal() = vec;
      v.meanCurvature(curvature);
      return curvature;
    }
    
    Halfedge_around_vertex_circulator hec = v.vertex_begin();
    Halfedge_around_vertex_circulator d = hec;
    if(hec == NULL)
    {
      std::cerr << "** isolated vertex" << std::endl;
      return 0.0; // isolated vertex
    }
      
    int degree = 0;
    CGAL_For_all(hec, d)
    {
      Vector ap = (hec->vertex()->point()-
                   hec->next()->next()->vertex()->point());
      double cot_beta  = cotangent(hec->next());
      double cot_alpha = cotangent(hec->opposite()->next());
      vec = vec +  (cot_alpha+cot_beta) * ap;
      degree++;
    }

    // mixed area
    double mixed = mixed_area(v,pure_voronoi);
    if(mixed > 0.0)
      vec = vec/mixed;
    else
      vec = CGAL::NULL_VECTOR;
    // CGAL_assertion(mixed > 0.0);

    // take norm
    double curvature = std::sqrt(vec*vec);
    v.meanCurvature(curvature);
    if(curvature != 0.0)
      v.meanCurvatureNormal() = vec/curvature;
    else
      v.meanCurvatureNormal() = vec;

    return curvature;
  }

  // return cotan of corner specified by he->vertex()
  double cotangent(Halfedge_handle he)
  {
    Vector pa = (he->next()->next()->vertex()->point()-he->vertex()->point());
    Vector pb = (he->next()->vertex()->point()-he->vertex()->point());
    // (pa . pb)/((pa x pb).len)
    double dot = (pa*pb);
    Vector cross = CGAL::cross_product(pa,pb);
    double mag = std::sqrt(cross*cross);
    if(mag != 0.0)
      return (dot/mag);
    else
      return 0.0; // undefined
  }

  // area mixed
  // see http://www.multires.caltech.edu/pubs/diffGeoOps.pdf, page 10
  double mixed_area(Vertex &v,
                    bool& pure_voronoi)
  {
    // voronoi area
    double voronoiArea = 0.0;
    Halfedge_around_vertex_circulator hec = v.vertex_begin();
    Halfedge_around_vertex_circulator d = hec;
    CGAL_assertion(hec != NULL);

    // for each edge
    CGAL_For_all(hec, d)
    {
      Vector ap = (hec->vertex()->point()-hec->next()->next()->vertex()->point());
      Vector pa = (hec->next()->next()->vertex()->point()-hec->vertex()->point());
      Vector pb = (hec->next()->vertex()->point()-hec->vertex()->point());
      Vector ab = (hec->next()->vertex()->point()-hec->next()->next()->vertex()->point());

      // one obtuse triangle, switch to other measure
      if((pa*pb) < 0.0 ||
         (ap*ab) < 0.0 ||
         (pb*ab) < 0.0)
      {
        pure_voronoi = false;
        return mixed_area_obtuse(v);
      }

      // continue
      double cot_beta  = cotangent(hec->next());
      double cot_alpha = cotangent(hec->opposite()->next());
      CGAL_assertion(hec->next() != hec->opposite()->next());
      voronoiArea += (cot_alpha+cot_beta) * (ap*ap);
    }
    voronoiArea *= 0.125; // 1/8*sum(cot+cot)(xi-xj)^2
    pure_voronoi = true;
    return voronoiArea;
  }

  // area mixed (obtuse case)
  // see http://www.multires.caltech.edu/pubs/diffGeoOps.pdf, page 10
  double mixed_area_obtuse(Vertex &v)
  {
    // voronoi area
    double area = 0.0;
    Halfedge_around_vertex_circulator hec = v.vertex_begin();
    Halfedge_around_vertex_circulator d = hec;

    // for each edge
    CGAL_For_all(hec, d)
    {
      Vector ap = (hec->vertex()->point()-hec->next()->next()->vertex()->point());
      Vector pa = (hec->next()->next()->vertex()->point()-hec->vertex()->point());
      Vector pb = (hec->next()->vertex()->point()-hec->vertex()->point());
      Vector ab = (hec->next()->vertex()->point()-hec->next()->next()->vertex()->point());

      // one obtuse triangle
      if((pa*pb) < 0.0 ||
         (ap*ab) < 0.0 ||
         (pb*ab) < 0.0)
      {
        Vector cross = CGAL::cross_product(pa,pb);
        double len = std::sqrt(cross*cross);
        area += (len/6.0);
      }
      else
      {
        // continue
        double cot_beta  = cotangent(hec->next());
        double cot_alpha = cotangent(hec->next()->next());
        area += 0.125*(cot_alpha * (pb*pb) +
                       cot_beta  * (ap*ap));
      }
    }
    return area;
  }

  // get any sharp or seaming halfedge with tag
  Halfedge_handle get_sharp_or_seaming_halfedge_tag(int tag)
  {
    // high priority is given to seaming halfedges
    Halfedge_iterator pHalfedge;
    for(pHalfedge = halfedges_begin();
        pHalfedge != halfedges_end();
        pHalfedge++)
    {
      if(pHalfedge->is_seaming() &&
         pHalfedge->tag() == tag)
        return pHalfedge;
    }
    // lower priority is given to sharp halfedges
    for(pHalfedge = halfedges_begin();
        pHalfedge != halfedges_end();
        pHalfedge++)
    {
      if(pHalfedge->is_sharp() &&
         pHalfedge->tag() == tag)
        return pHalfedge;
    }
    return NULL;
  }

  // get any border halfedge with tag
  Halfedge_handle get_border_halfedge_tag(int tag)
  {
    Halfedge_iterator pHalfedge;
    for(pHalfedge = halfedges_begin();
        pHalfedge != halfedges_end();
        pHalfedge++)
    {
      if(pHalfedge->is_border() &&
         pHalfedge->tag() == tag)
        return pHalfedge;
    }
    return NULL;
  }

  // get any sharp or seaming halfedge attached to a corner with tag
  Halfedge_handle get_border_halfedge_tag(Vertex_handle pVertex,
                                          int tag)
  {
    Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator d = pHalfEdge;
    CGAL_For_all(pHalfEdge,d)
    {
      Halfedge_handle next = pHalfEdge->next();
      if(next->is_border() && next->tag() == tag)
        return next;
    }
    return NULL;
  }

  // get any sharp or seaming halfedge attached to a corner
  Halfedge_handle get_border_halfedge(Vertex_handle pVertex)
  {
    Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator d = pHalfEdge;
    CGAL_For_all(pHalfEdge,d)
      if(pHalfEdge->next()->is_border())
        return pHalfEdge->next();
    return NULL;
  }

  // get any sharp or seaming halfedge attached to a corner with tag
  Halfedge_handle get_sharp_or_seaming_halfedge_tag(Vertex_handle pVertex,
                                                    int tag)
  {
    Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator d = pHalfEdge;
    CGAL_For_all(pHalfEdge,d)
      if((pHalfEdge->is_sharp() || pHalfEdge->is_seaming()) &&
         pHalfEdge->tag() == tag)
        return pHalfEdge;
    return NULL;
  }

  // get any corner with a sharp/seaming halfedge with tag
  Vertex_handle get_corner_sharp_or_seaming_halfedge_tag(int tag)
  {
    Vertex_iterator pVertex;
    for(pVertex = vertices_begin();
        pVertex != vertices_end();
        pVertex++)
    {
      if(pVertex->is_corner() &&
         (get_sharp_or_seaming_halfedge_tag(pVertex,tag)) != NULL)
        return pVertex;
    }
    return NULL;
  }

  // get any corner with a border halfedge with tag
  Vertex_handle get_corner_border_halfedge_tag(int tag)
  {
    Vertex_iterator pVertex;
    for(pVertex = vertices_begin();
        pVertex != vertices_end();
        pVertex++)
    {
      if(pVertex->is_corner() &&
         (get_border_halfedge_tag(pVertex,tag)) != NULL)
        return pVertex;
    }
    return NULL;
  }


  // get index of the longest backbone
  int get_index_longest_backbone()
  {
    int index = 0;
    double max = 0.0;
    // #backbones
    int nb = (*m_skeleton.backbones()).size();
    for(int i=0;i<nb;i++)
    {
      backbone *pBackbone = (*m_skeleton.backbones())[i];
      double length = len(pBackbone);
      if(length>max)
      {
        index = i;
        max = length;
      }
    }
    return index;
  }

  // count #boundaries
  // return the number of boundary backbones
  int nb_boundaries()
  {
    int nb = 0;
    tag_halfedges(0);
    Halfedge_handle seed_halfedge = NULL;
    while((seed_halfedge = get_border_halfedge_tag(0)) != NULL)
    {
      nb++;
      seed_halfedge->tag(1);
      Vertex_handle seed_vertex = seed_halfedge->prev()->vertex();
      Halfedge_handle current_halfedge = seed_halfedge;
      Halfedge_handle next_halfedge;
      do
      {
        next_halfedge = current_halfedge->next();
        next_halfedge->tag(1);
        current_halfedge = next_halfedge;
      }
      while(next_halfedge->prev()->vertex() != seed_vertex);
    }
    return nb;
  }

  // compute the genus
  // V - E + F + B = 2 (C - G)
  // C -> #connected components
  // G : genus
  // B : #boundaries
  int genus()
  {
    int b = nb_boundaries();
    int v = size_of_vertices();
    int e = size_of_halfedges()/2;
    int f = size_of_facets();
    int genus = (2+e-b-f-v)/2;
    std::cerr << "  " << v << " vertices, " << f << " faces, ";
    std::cerr << e << " edges, " << b << " boundary(ies), genus " << genus << std::endl;
    return genus;
  }

  // get any adjacent halfedge with target diff than pVertexDiff
  Halfedge_handle get_halfedge_target_diff(Vertex_handle pVertex,
                                           Vertex_handle pVertexDiff)
  {
    Polyhedron_ex::Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
    CGAL_assertion(pHalfEdge != NULL);
    Polyhedron_ex::Halfedge_around_vertex_circulator d = pHalfEdge;
    CGAL_For_all(pHalfEdge,d)
      if(pHalfEdge->opposite()->vertex() != pVertexDiff)
        return pHalfEdge;
    return NULL;
  }

  Vertex_handle search_farthest(Vertex_handle vertex)
  {
    // pick one arbitrary seed vertex
    double max = 0.0;
    Vertex_handle farthest = NULL;
    Vertex_iterator v;
    for(v = vertices_begin();
        v != vertices_end();
        v++)
    {
      double distance = My_kernel::len(v->point()-vertex->point());
      if(distance > max)
      {
        max = distance;
        farthest = v;
      }
    }
    CGAL_assertion(farthest != NULL);
    return farthest;
  }

  // close one arbitrary hole
  // return true if closed
  bool close_one_hole()
  {
    Halfedge_iterator pHalfedge;
    for(pHalfedge = halfedges_begin();
        pHalfedge != halfedges_end();
        pHalfedge++)
    {
      if(pHalfedge->is_border())
      {
        Halfedge_handle he = pHalfedge;
        fill_hole(he);
        std::cerr << "  close hole" << std::endl;
        return true;
      }
    }
    return false;
  }

  // close all holes
  void close_holes()
  {
    while(close_one_hole());
  }

  // compute  total len of a backbone
  double len(backbone *pBackbone)
  {
    std::list<Polyhedron_ex::Halfedge_handle> *pHalfedges = pBackbone->halfedges();
    std::list<Polyhedron_ex::Halfedge_handle>::iterator pHalfedge;
    double len = 0.0;
    for(pHalfedge = pHalfedges->begin();
        pHalfedge != pHalfedges->end();
        pHalfedge++)
     {
       Polyhedron_ex::Halfedge_handle he = (*pHalfedge);
       Vector v = (he->vertex()->point()-he->prev()->vertex()->point());
       len += My_kernel::len(v);
     }
     return len;
   }

   // compute distance from facet center to halfedge center
   double distance(Facet_handle pFacet,
                   Halfedge_handle pHalfedge)
   {
     // we assume
     Point center_facet = pFacet->center();

    Vector v = (pHalfedge->opposite()->vertex()->point()
                -pHalfedge->vertex()->point());
     Point center_halfedge = pHalfedge->vertex()->point() + (v/2);
    Vector d = center_facet-center_halfedge;
    return My_kernel::len(d);
   }

  // set is_inner flag to all vertices from a backbone
  void set_is_inner_param_vertices(backbone *pBackbone,
                                   bool is_inner)
  {
    assert(pBackbone != NULL);
    std::list<Polyhedron_ex::Halfedge_handle> *pHalfedges = pBackbone->halfedges();
    assert(pHalfedges != NULL);
    assert(pHalfedges->size() > 0);
    std::list<Polyhedron_ex::Halfedge_handle>::iterator pHalfedge;
    for(pHalfedge = pHalfedges->begin();
        pHalfedge != pHalfedges->end();
        pHalfedge++)
     {
       Polyhedron_ex::Halfedge_handle he = (*pHalfedge);
      assert(he != NULL);
       he->vertex()->is_inner_param(is_inner);
       he->opposite()->vertex()->is_inner_param(is_inner);
     }
   }

   void farthest_point_aligned(Vertex_handle &pVertexMin,
                               Vertex_handle &pVertexMax)
   {
     double xmin,xmax,ymin,ymax,zmin,zmax;
    Vertex_handle pVertex_xMin = vertex_min(0,xmin);
    Vertex_handle pVertex_xMax = vertex_max(0,xmax);
    Vertex_handle pVertex_yMin = vertex_min(1,ymin);
    Vertex_handle pVertex_yMax = vertex_max(1,ymax);
    Vertex_handle pVertex_zMin = vertex_min(2,zmin);
    Vertex_handle pVertex_zMax = vertex_max(2,zmax);
    double xdiff = xmax-xmin;
    double ydiff = ymax-ymin;
    double zdiff = zmax-zmin;
    double max = std::max(std::max(xdiff,ydiff),zdiff);
    if(max == xdiff)
    {
      pVertexMin = pVertex_xMin;
      pVertexMax = pVertex_xMax;
    }
    else
      if(max == ydiff)
      {
        pVertexMin = pVertex_yMin;
        pVertexMax = pVertex_yMax;
      }
      else
      {
        pVertexMin = pVertex_zMin;
        pVertexMax = pVertex_zMax;
      }
   }

   // deduce uv from barycentric coordinates
   void uv_from_barycentric_coordinates(Facet_handle &pFacet,
                                        double &u,
                                        double &v,
                                        const double &alpha,
                                        const double &beta,
                                        const double &gamma)
   {
    Halfedge_around_facet_circulator h = pFacet->facet_begin();
    double u0 = h->u();
    double v0 = h->v();
    double u1 = h->next()->u();
    double v1 = h->next()->v();
    double u2 = h->next()->next()->u();
    double v2 = h->next()->next()->v();
    u = alpha*u0 + beta*u1 + gamma*u2;
    v = alpha*v0 + beta*v1 + gamma*v2;
  }

   // barycentric coordinates in a triangle
   // return true if inside
   void barycentric_coordinates(double u,
                                double v,
                                Facet_handle &pFacet,
                                double &alpha,
                                double &beta,
                                double &gamma,
                                bool fix)
   {
    Halfedge_around_facet_circulator h = pFacet->facet_begin();
    double u0 = h->u();
    double v0 = h->v();
    double u1 = h->next()->u();
    double v1 = h->next()->v();
    double u2 = h->next()->next()->u();
    double v2 = h->next()->next()->v();

    // total triangle area
    double a = (u1-u0)*(v2-v0)-(v1-v0)*(u2-u0);
    if(a == 0)
    {
      alpha = 1;
      beta  = 0;
      gamma = 0;
      std::cerr << "null area" << std::endl;
      return;
    }

    double inva = 1.0/a;
    double a0 = (u1-u)*(v2-v)-(v1-v)*(u2-u);
    double a1 = (u2-u)*(v0-v)-(v2-v)*(u0-u);
    alpha = inva*a0;
    beta  = inva*a1;
    gamma = 1.0-alpha-beta;

    // optional fix
    if(fix)
    {
      //TRACE("** fix (%g;%g;%g) -> ",alpha,beta,gamma);
      alpha = (alpha < 0) ? 0 : (alpha > 1) ? 1 : alpha;
      beta  = (beta  < 0) ? 0 : (beta  > 1) ? 1 : beta;
      gamma = 1.0-alpha-beta;
      if(gamma < 0) // BIG robustness issue
      {
        gamma = 0;
        beta  = 1-alpha;
      }
      //TRACE("(%g;%g;%g)\n",alpha,beta,gamma);
    }
   }

   // return true if uv coordinates are
   // inside the facet
   bool is_inside(double u,
                  double v,
                  Facet_handle &pFacet,
                  double &alpha,
                  double &beta,
                  double &gamma)
   {
     // compute barycentric coordinates
     barycentric_coordinates(u,v,pFacet,alpha,beta,gamma,false);

    // inside
    if(alpha >= 0 &&
       alpha <= 1 &&
       beta  >= 0 &&
       beta  <= 1 &&
       gamma >= 0 &&
       gamma <= 1)
      return true;
    else
      return false;
  }

  // locate facet
  Facet_handle locate_facet(double u,
                            double v,
                            double &alpha,
                            double &beta,
                            double &gamma,
                            Facet_handle &pFacetTryFirst)
  {
    // try this one first
    if(pFacetTryFirst != NULL)
      if(is_inside(u,v,pFacetTryFirst,alpha,beta,gamma))
        return pFacetTryFirst;

     Facet_iterator pFacet;
     for(pFacet = facets_begin();
         pFacet != facets_end();
         pFacet++)
     {
      if(is_inside(u,v,pFacet,alpha,beta,gamma))
        return pFacet;
    }
    return NULL;
  }

    // get any inner vertex
    Vertex_handle get_any_inner_vertex()
    {
      Vertex_iterator pVertex;
      for(pVertex = vertices_begin();
          pVertex != vertices_end();
          pVertex++)
        if(!is_border(*pVertex))
          return pVertex;
     return NULL;
    }
    // get any halfedge with tag diff
    Halfedge_handle get_any_halfedge_tag_diff(int tag)
    {
      Halfedge_iterator pHalfedge;
      for(pHalfedge = halfedges_begin();
          pHalfedge != halfedges_end();
          pHalfedge++)
        if(pHalfedge->tag() != tag)
          return pHalfedge;
     return NULL;
    }
    // get any halfedge with tag
    Halfedge_handle get_any_halfedge_tag(int tag)
    {
      Halfedge_iterator pHalfedge;
      for(pHalfedge = halfedges_begin();
          pHalfedge != halfedges_end();
          pHalfedge++)
        if(pHalfedge->tag() == tag)
          return pHalfedge;
     return NULL;
    }
}; // end class PolyhedronEx



#endif // _CGAL_DEFS_



