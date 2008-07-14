#ifndef _MAKE_QUAD_
#define _MAKE_QUAD_

#include <CGAL/Polyhedron_incremental_builder_3.h>

template <class HDS,class Polyhedron,class Kernel>
class CModifierQuad : public CGAL::Modifier_base<HDS>
{
private:
	typedef typename Kernel::Point_3 Point;
	typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> builder;
	Point m_points[4];

public:

	// life cycle
	CModifierQuad(const Point& a,
		      const Point& b,
		      const Point& c,
		      const Point& d)
	{
		m_points[0] = a;
		m_points[1] = b;
		m_points[2] = c;
		m_points[3] = d;
	}
	~CModifierQuad() {}

	void operator()( HDS& hds)
	{
		builder B(hds,true);
		B.begin_surface(4,1,8);
		
		for(int i=0;i<4;i++)
			B.add_vertex(m_points[i]);

		B.begin_facet();
			B.add_vertex_to_facet(0);
			B.add_vertex_to_facet(1);
			B.add_vertex_to_facet(2);
			B.add_vertex_to_facet(3);
		B.end_facet();
		
		B.end_surface();
	}
};

template <class Polyhedron,class Kernel>
class Make_quad
{
public:
	typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
	typedef typename Kernel::Point_3 Point;
	Make_quad() {}
	~Make_quad() {}

public:
	void run(const Point& a,
		 const Point& b,
		 const Point& c,
		 const Point& d,
		 Polyhedron &output)
	{
		CModifierQuad<HalfedgeDS,Polyhedron,Kernel> quad(a,b,c,d);
		output.delegate(quad);
	}
};

#endif // _MAKE_QUAD_
