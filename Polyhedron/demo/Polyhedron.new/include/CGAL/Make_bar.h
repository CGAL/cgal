#ifndef _MAKE_BAR_
#define _MAKE_BAR_

#include <CGAL/Polyhedron_incremental_builder_3.h>

template <class HDS,class Polyhedron,class Kernel>
class CModifierBar : public CGAL::Modifier_base<HDS>
{
private:
	typedef typename Kernel::Point_3 Point;
	typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> builder;
	Point m_points[8];

public:

	// life cycle
	CModifierBar(Point points[8])
	{
		for(int i=0;i<8;i++)
			m_points[i] = points[i];
	}
	~CModifierBar() {}

	void operator()( HDS& hds)
	{
		builder B(hds,true);
		B.begin_surface(3,1,6);
		
		for(int i=0;i<8;i++)
			B.add_vertex(m_points[i]);

		B.begin_facet();
			B.add_vertex_to_facet(0);
			B.add_vertex_to_facet(1);
			B.add_vertex_to_facet(2);
			B.add_vertex_to_facet(3);
		B.end_facet();

		B.begin_facet();
			B.add_vertex_to_facet(0);
			B.add_vertex_to_facet(4);
			B.add_vertex_to_facet(5);
			B.add_vertex_to_facet(1);
		B.end_facet();

		B.begin_facet();
			B.add_vertex_to_facet(1);
			B.add_vertex_to_facet(5);
			B.add_vertex_to_facet(6);
			B.add_vertex_to_facet(2);
		B.end_facet();

		B.begin_facet();
			B.add_vertex_to_facet(2);
			B.add_vertex_to_facet(6);
			B.add_vertex_to_facet(7);
			B.add_vertex_to_facet(3);
		B.end_facet();

		B.begin_facet();
			B.add_vertex_to_facet(3);
			B.add_vertex_to_facet(7);
			B.add_vertex_to_facet(4);
			B.add_vertex_to_facet(0);
		B.end_facet();

		B.begin_facet();
			B.add_vertex_to_facet(4);
			B.add_vertex_to_facet(7);
			B.add_vertex_to_facet(6);
			B.add_vertex_to_facet(5);
		B.end_facet();

		B.end_surface();
	}
};

template <class Polyhedron,class Kernel>
class Make_bar
{
public:
	typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
	typedef typename Kernel::Point_3 Point;
	Make_bar() {}
	~Make_bar() {}

public:
	void run(Point points[8],
		 Polyhedron &output)
	{
		CModifierBar<HalfedgeDS,Polyhedron,Kernel> bar(points);
		output.delegate(bar);
	}
};

#endif // _MAKE_BAR_
