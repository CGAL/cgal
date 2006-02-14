///////////////////////////////////////////////////////////////////////////
//																																			 //
//	Class: Enriched_polyhedron                                           //
//																																			 //
///////////////////////////////////////////////////////////////////////////

#ifndef	_POLYGON_MESH_
#define	_POLYGON_MESH_

// CGAL	stuff
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <list>

// a refined facet with a normal 
template <class	Refs,	class	T, class P,	class	Norm>
class	Enriched_facet : public	CGAL::HalfedgeDS_face_base<Refs, T>
{
	// normal
	Norm m_normal;

public:

	// life	cycle
	Enriched_facet()
	{
	}

	// normal
	typedef	Norm Normal_3;
	Normal_3&	normal() { return	m_normal;	}
	const	Normal_3&	normal() const { return	m_normal;	}
};

// a refined halfedge with a general tag 
template <class	Refs,	class	Tprev, class Tvertex,	class	Tface, class Norm>
class	Enriched_halfedge	:	public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
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
};

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


// Enriched polyhedron
template <class	kernel,	class	items>
class	Enriched_polyhedron	:	public CGAL::Polyhedron_3<kernel,items>
{
public :
	typedef	typename kernel::FT	FT;
	typedef	typename kernel::Point_3 Point;

public :

	// life	cycle
	Enriched_polyhedron()	
	{
	}
	virtual	~Enriched_polyhedron() 
	{
	}

	// tag all halfedges
	void tag_halfedges(const int tag)
	{
		for(Halfedge_iterator pHalfedge	=	halfedges_begin();
				pHalfedge	!= halfedges_end();
				pHalfedge++)
			pHalfedge->tag(tag);
	}

	// draw	edges
	void gl_draw_edges()
	{
		::glBegin(GL_LINES);

		Halfedge_iterator it;
		for(it = edges_begin();
			  it != edges_end();
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

	// draw	sites
	void gl_draw_sites(const unsigned char r,
		                 const unsigned char g,
										 const unsigned char b,
										 const float size)
	{
		::glPointSize(size);
		::glColor3ub(r,g,b);
		::glBegin(GL_POINTS);
		Point_iterator it;
		for(it = points_begin();
			  it != points_end();
				it++)
		{
			const Point& p = *it;
			::glVertex3d(p[0],p[1],p[2]);
		}
		::glEnd();
	}

	// normals (per	facet, then	per	vertex)
	void compute_normals_per_facet()
	{
		std::for_each(facets_begin(),facets_end(),Facet_normal());
	}
	void compute_normals_per_vertex()
	{
		std::for_each(vertices_begin(),vertices_end(),Vertex_normal());
	}
	void compute_normals()
	{
		compute_normals_per_facet();
		compute_normals_per_vertex();
	}
};

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
			typename Facet::Normal_3 normal	=	CGAL::cross_product(
				h->next()->vertex()->point() - h->vertex()->point(),
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
			TRACE("degenerate face\n");
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
				Vertex::Halfedge_around_vertex_const_circulator	pHalfedge	=	v.vertex_begin();
				Vertex::Halfedge_around_vertex_const_circulator	begin	=	pHalfedge;
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

#endif //	_POLYGON_MESH_
