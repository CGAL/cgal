///////////////////////////////////////////////////////////////////////////
//
// Class: Enriched_polyhedron
//
//
///////////////////////////////////////////////////////////////////////////

#ifndef	CGAL_POLYGON_MESH_
#define	CGAL_POLYGON_MESH_

// CGAL	stuff
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <list>

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// compute facet normal
struct Facet_normal	// (functor)
{
  template <class	Facet>
    void operator()(Facet& f)
    {
      typename Facet::Normal_3 sum = CGAL::NULL_VECTOR;
      typename Facet::Halfedge_around_facet_circulator h = f.facet_begin();
      do
	{
	  typename Facet::Normal_3 normal =
	    CGAL::cross_product( h->next()->vertex()->point() - h->vertex()->point(),
				 h->next()->next()->vertex()->point() - h->next()->vertex()->point());
	  double sqnorm	=	normal * normal;
	  if(sqnorm	!= 0)
	    normal = normal	/	(float)std::sqrt(sqnorm);
	  sum	=	sum	+	normal;
	}
      while(++h	!= f.facet_begin());
      float	sqnorm = sum * sum;
      if(sqnorm	!= 0.0)
	f.normal() = sum / std::sqrt(sqnorm);
      else
	{
	  f.normal() = CGAL::NULL_VECTOR;
	}
    }
};


// compute vertex	normal
struct Vertex_normal //	(functor)
{
  template <class	Vertex>
    void operator()(Vertex&	v)
    {
      typename Vertex::Normal_3	normal = CGAL::NULL_VECTOR;
      typename Vertex::Halfedge_around_vertex_const_circulator
	pHalfedge = v.vertex_begin();
      typename Vertex::Halfedge_around_vertex_const_circulator
	begin =	pHalfedge;

      CGAL_For_all(pHalfedge,begin)
	if(!pHalfedge->is_border())
	  normal = normal	+	pHalfedge->facet()->normal();
      float	sqnorm = normal * normal;
      if(sqnorm != 0.0f)
	v.normal() = normal	/	(float)std::sqrt(sqnorm);
      else
	v.normal() = CGAL::NULL_VECTOR;
    }
};

// a refined facet with a normal
template <class	Refs,	class	T, class P,	class	Norm>
class	Enriched_facet : public	CGAL::HalfedgeDS_face_base<Refs, T>
{
  // normal
  Norm  m_normal;
  int	m_tag;

public:

  // life	cycle
  Enriched_facet()
    {
      m_tag = 0;
    }

  // normal
  typedef	Norm Normal_3;
  Normal_3&	normal() { return	m_normal;	}
  const	Normal_3&	normal() const { return	m_normal;	}

  //tag
  int& tag() { return m_tag; }
  const int& tag() const { return m_tag; }
  void tag(const int& t) { m_tag = t; }
};

//-----------------------------------------------------------------------------
// a refined halfedge with a general tag
template <class	Refs,	class	Tprev, class Tvertex,	class	Tface, class Norm>
class	Enriched_halfedge	:	public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
public :
  float r,v,b;
private:

  // general-purpose tag
  int m_tag;

public:

  // life	cycle
  Enriched_halfedge()
    {
    }

  // tag
  int& tag() { return m_tag; }
  const int& tag() const { return m_tag; }
  void tag(const int& t) { m_tag = t; }

};


//-----------------------------------------------------------------------------
// a refined vertex with a normal and a tag
template <class	Refs,	class	T, class P,	class	Norm>
class	Enriched_vertex	:	public CGAL::HalfedgeDS_vertex_base<Refs,	T, P>
{
  // normal
  Norm m_normal;

public:
  // life	cycle
  Enriched_vertex()	 {}
  // repeat	mandatory	constructors
  Enriched_vertex(const	P& pt)
    :	CGAL::HalfedgeDS_vertex_base<Refs, T,	P>(pt)
    {
    }

  // normal
  typedef	Norm Normal_3;
  Normal_3&	normal() { return	m_normal;	}
  const	Normal_3&	normal() const { return	m_normal;	}

  // number of vertex
  unsigned int num;
};

//-----------------------------------------------------------------------------
// A redefined items class for the Polyhedron_3
// with	a	refined	vertex class that	contains a
// member	for	the	normal vector	and	a	refined
// facet with	a	normal vector	instead	of the
// plane equation	(this	is an	alternative
// solution	instead	of using
// Polyhedron_traits_with_normals_3).

struct Enriched_items	:	public CGAL::Polyhedron_items_3
{
  // wrap	vertex
  template <class	Refs,	class	Traits>
  struct Vertex_wrapper
  {
    typedef	typename Traits::Point_3	Point;
    typedef	typename Traits::Vector_3	Normal;
    typedef	Enriched_vertex<Refs,
      CGAL::Tag_true,
      Point,
      Normal>	Vertex;
  };

  // wrap	face
  template <class	Refs,	class	Traits>
  struct Face_wrapper
  {
    typedef	typename Traits::Point_3	Point;
    typedef	typename Traits::Vector_3	Normal;
    typedef	Enriched_facet<Refs,
      CGAL::Tag_true,
      Point,
      Normal> Face;
  };

  // wrap	halfedge
  template <class	Refs,	class	Traits>
  struct Halfedge_wrapper
  {
    typedef	typename Traits::Vector_3	Normal;
    typedef	Enriched_halfedge<Refs,
      CGAL::Tag_true,
      CGAL::Tag_true,
      CGAL::Tag_true,
      Normal>	Halfedge;
  };
};

//-----------------------------------------------------------------------------
// Enriched polyhedron
template <class	kernel,	class	items>
class	Enriched_polyhedron:public CGAL::Polyhedron_3<kernel,items>
{
public :
  typedef typename kernel::FT	FT;
  typedef typename kernel::Point_3 Point;
  typedef typename kernel::Vector_3 Vector;
  typedef typename Enriched_polyhedron<kernel, items>::Vertex_handle Vertex_handle;
  typedef typename Enriched_polyhedron<kernel, items>::Facet_handle Facet_handle;

  typedef typename Enriched_polyhedron<kernel, items>::Facet Facet;

  typedef typename Enriched_polyhedron<kernel, items>::Face_handle Face_handle;
  typedef typename Enriched_polyhedron<kernel, items>::Halfedge_handle Halfedge_handle;
  typedef typename Enriched_polyhedron<kernel, items>::Facet_iterator Facet_iterator;
  typedef typename Enriched_polyhedron<kernel, items>::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
  typedef typename Enriched_polyhedron<kernel, items>::Halfedge_iterator Halfedge_iterator;
  typedef typename Enriched_polyhedron<kernel, items>::Point_iterator Point_iterator;
  typedef typename Enriched_polyhedron<kernel, items>::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;


typedef typename Enriched_polyhedron<kernel, items>::Vertex_iterator Vertex_iterator;
typedef typename Enriched_polyhedron<kernel, items>::Vertex Vertex;

typedef typename Enriched_polyhedron<kernel, items>::Edge_iterator Edge_iterator ;
public :

  // life	cycle
  Enriched_polyhedron()
    {
    }
  virtual	~Enriched_polyhedron()
    {
    }



  bool is_border(Vertex_handle v)
    {
      Halfedge_around_vertex_circulator	he	=	v->vertex_begin();
      if(he ==	NULL)	// isolated	vertex
	return true;

      Halfedge_around_vertex_circulator	cir	=	he;
      CGAL_For_all(cir,he)
	{
	  if(cir->is_border())
	    return true;
	}
      return false;
    }

  // draw	edges
  void gl_draw_edges()
    {
      ::glColor3f(1,0,0);
      ::glBegin(GL_LINES);

      Halfedge_iterator it;
      for(it = this->edges_begin();
	  it != this->edges_end();
	  it++)
	{
	  Halfedge_handle he = it;
	  const Point& p1 =	he->opposite()->vertex()->point();
	  const Point& p2 =	he->vertex()->point();
	  ::glVertex3d(p1[0],p1[1],p1[2]);
	  ::glVertex3d(p2[0],p2[1],p2[2]);
	}
      ::glEnd();
     }


  // draw vertices
  void gl_draw_vertices()
    {
      ::glPointSize(5.0f);
      ::glBegin(GL_POINTS);
      Point_iterator it;
      for(it = this->points_begin();
	  it != this->points_end();
	  it++)
	{
	  const Point& p = *it;
	  ::glVertex3d(p[0],p[1],p[2]);
	}
      ::glEnd();
    }

  // draw facets
  void gl_draw_facets(const bool smooth)
    {
      glEnable(GL_LIGHTING);
      static GLfloat agray[4] = {1,1,1, 1.0 };
      glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, agray);

     glEnable(GL_POLYGON_OFFSET_FILL);
     glPolygonOffset( 1.0, 1.0 );
    glPolygonMode(GL_FRONT, GL_FILL);

      Facet_iterator hFacet;
      for(hFacet = this->facets_begin();
	  hFacet != this->facets_end();
	  hFacet++)
	gl_draw_facet(hFacet,smooth);

      glPolygonOffset( 0.0, 0.0 );
      glDisable(GL_POLYGON_OFFSET_FILL);

   }

  void gl_draw_facet(Facet_handle hFacet,
		     const bool smooth)
    {

      ::glBegin(GL_POLYGON);

      // one normal per facet
      if(!smooth)
	{
	  const typename Facet::Normal_3& normal 	= hFacet->normal();
	  ::glNormal3f(normal[0],normal[1],normal[2]);
	}

      Halfedge_around_facet_circulator he = hFacet->facet_begin();
      do
	{
	  //one normal per vertices
	  if(smooth)
	    {
	      const typename Vertex::Normal_3& normal	=	he->vertex()->normal();
	      ::glNormal3f(normal[0],normal[1],normal[2]);
	    }
	  const Point& p =	he->vertex()->point();
	  ::glVertex3d(p[0],p[1],p[2]);
	}
      while(++he	!= hFacet->facet_begin());
      ::glEnd();
    }


  Point iso_barycentre ( Face_handle hFacet )
    {
      // calcul de l'isobarycentre de la face
      double x,y,z;
      x=0;
      y=0;
      z=0;
      int n = 0;

      Halfedge_around_facet_circulator he = hFacet->facet_begin();
      do
	{
	  const Point& p =	he->vertex()->point();
	  x = x + p[0];
	  y = y + p[1];
	  z = z + p[2];
	  n++;
	}
      while(++he	!= hFacet->facet_begin());

      return Point( FT(x/n), FT(y/n),FT(z/n));
    }

  void gl_draw_facets_normal()
    {
      ::glColor3f(1,0,0);
      ::glBegin(GL_LINES);

      Facet_iterator hFacet;
      for(hFacet = this->facets_begin();
	  hFacet != this->facets_end();
	  hFacet++)
	{
	  Point p = iso_barycentre(hFacet);
	  const typename Facet::Normal_3& normal = hFacet->normal();
	  typename kernel::Vector_3 n = normal /
	    sqrt( pow(normal[0],2)+pow(normal[1],2)+pow(normal[2],2));
	  ::glVertex3d(p[0],p[1],p[2]);
	  ::glVertex3d(p[0]+n[0],p[1]+n[1],p[2]+n[2]);

	}

      ::glEnd();
    }


  void gl_draw_vertices_normal()
    {
      ::glBegin(GL_LINES);
      Vertex_iterator it;
      for(it = this->vertices_begin();
	  it != this->vertices_end();
	  it++)
	{
	  const typename Vertex::Normal_3& normal	=	it->normal();
	  const Point& p = it->point();
	  typename kernel::Vector_3 n = normal /
	    sqrt( pow(normal[0],2)+pow(normal[1],2)+pow(normal[2],2));
	  ::glVertex3d(p[0],p[1],p[2]);
	  ::glVertex3d(p[0]+n[0],p[1]+n[1],p[2]+n[2]);
	}
      ::glEnd();
    }

  // normals (per	facet, then	per	vertex)
  void compute_normals_per_facet()
    {
      std::for_each(this->facets_begin(),this->facets_end(),Facet_normal());
    }
void compute_normals_per_vertex()
{
  std::for_each(this->vertices_begin(),this->vertices_end(),Vertex_normal());
}
void compute_normals()
{
  compute_normals_per_facet();
  compute_normals_per_vertex();
}


// count number of componant
int nb_component()
{
  unsigned int nb = 0;
  tag_facets(0);
  Facet_handle seed_facet	=	NULL;
  while((seed_facet	=	this->get_facet_tag(0))	!= NULL)
    {
      nb++;
      tag_component(seed_facet,0,1);
    }
  return nb;
}

// tag all facets
void tag_facets(const	int	tag)
{
  for(Facet_iterator pFace	=	this->facets_begin();
      pFace	!= this->facets_end();
      pFace++)
    pFace->tag(tag);
}

// tag component
void tag_component(Facet_handle	pSeedFacet,
		   const int tag_free,
		   const int tag_done)
{
  pSeedFacet->tag(tag_done);
  std::list<Facet_handle> facets;
  facets.push_front(pSeedFacet);
  while(!facets.empty())
    {
      Facet_handle pFacet = facets.front();
      facets.pop_front();
      pFacet->tag(tag_done);
      Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
      Halfedge_around_facet_circulator end = pHalfedge;
      CGAL_For_all(pHalfedge,end)
	{
	  Facet_handle pNFacet = pHalfedge->opposite()->facet();
	  if(pNFacet !=	NULL && pNFacet->tag() == tag_free)
	    {
	      facets.push_front(pNFacet);
	      pNFacet->tag(tag_done);
	    }
	}
    }
}


// tag all halfedges
void tag_halfedges(const int tag)
{
  for(Halfedge_iterator pHalfedge	=	this->halfedges_begin();
      pHalfedge	!= this->halfedges_end();
      pHalfedge++)
    pHalfedge->tag(tag);
}

// count #boundaries
unsigned int nb_boundaries()
{
  unsigned int nb = 0;
  tag_halfedges(0);
  Halfedge_handle	seed_halfedge	=	NULL;
  while((seed_halfedge = this->get_border_halfedge_tag(0)) !=	NULL)
    {
      nb++;
      seed_halfedge->tag(1);
      Vertex_handle	seed_vertex	=	seed_halfedge->prev()->vertex();
      Halfedge_handle	current_halfedge = seed_halfedge;
      Halfedge_handle	next_halfedge;
      do
	{
	  next_halfedge	=	current_halfedge->next();
	  next_halfedge->tag(1);
	  current_halfedge = next_halfedge;
	}
      while(next_halfedge->prev()->vertex()	!= seed_vertex);
    }
  return nb;
}

//------------------------------------------------------------------------------
void subdivide_sqrt3 ()
{
  // check for valid polygon mesh
  if(this->size_of_facets() == 0)
    return;


  // subdivision
  // We use that new vertices/halfedges/facets are appended at the end.
  std::size_t nv = this->size_of_vertices();
  Vertex_iterator last_v = this->vertices_end();
  -- last_v;  // the last of the old vertices
  Edge_iterator last_e = this->edges_end();
  -- last_e;  // the last of the old edges
  Facet_iterator last_f = this->facets_end();
  -- last_f;  // the last of the old facets

  Facet_iterator f = this->facets_begin();    // create new centre vertices
  do {
    star_facet( f);
  } while ( f++ != last_f);

  //Edge_iterator e = edges_begin();              // flip the old edges
  //++ last_e; // make it the past-the-end position again
  //while ( e != last_e) {
  //	Halfedge_handle h = e;
  //	++e; // careful, incr. before flip since flip destroys current edge
  //	flip_edge( h);
  //};
  CGAL_postcondition(this->is_valid());


}

void flip_edge(Halfedge_handle e)
{
  if(e->is_border_edge())
    return;
  Halfedge_handle h = e->next();
  join_facet( e);
  split_facet( h, h->next()->next());
}

void star_facet(Facet_iterator f)
{
  Vector vec( 0.0, 0.0, 0.0);
  std::size_t order = 0;
  Halfedge_around_facet_circulator h = f->facet_begin();
  do {
    vec = vec + ( h->vertex()->point() - CGAL::ORIGIN);
    ++ order;
  } while ( ++h != f->facet_begin());
  CGAL_assertion( order >= 3); // guaranteed by definition of Polyhedron
  Point centre =  CGAL::ORIGIN + (vec / (typename kernel::FT)order);

  Halfedge_handle new_centre = create_centre_vertex( f->halfedge());
  new_centre->vertex()->point() = centre;
}


















};



#endif //	CGAL_POLYGON_MESH_
