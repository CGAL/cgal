#ifndef _MAKE_QUAD_SOUP_
#define _MAKE_QUAD_SOUP_

#include <CGAL/Polyhedron_incremental_builder_3.h>

template <class HDS,class Polyhedron,class Kernel,class InputIterator>
class CModifierQuadSoup : public CGAL::Modifier_base<HDS>
{
private:
	typedef typename Kernel::Point_3 Point;
	typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> builder;
	InputIterator m_begin, m_end;

public:

	// life cycle
	CModifierQuadSoup(InputIterator begin,InputIterator end)
		: m_begin(begin), m_end(end)
	{
	}
	~CModifierQuadSoup() {}

	// make a quad soup
	void operator()( HDS& hds)
	{
		builder B(hds,true);
		B.begin_surface(4,1,8);

		int index = 0;
		InputIterator it;
		for(it = m_begin; it != m_end;)
		{
			B.add_vertex(*it++);
			B.add_vertex(*it++);
			B.add_vertex(*it++);
			B.add_vertex(*it++);

			B.begin_facet();
				B.add_vertex_to_facet(index++);
				B.add_vertex_to_facet(index++);
				B.add_vertex_to_facet(index++);
				B.add_vertex_to_facet(index++);
			B.end_facet();
		}
		B.end_surface();
	}
};

template <class Polyhedron,class Kernel,class InputIterator>
class Make_quad_soup
{
public:
	typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
	Make_quad_soup() {}
	~Make_quad_soup() {}

public:
	void run(InputIterator begin,
		 InputIterator end,
		 Polyhedron &output)
	{
		CModifierQuadSoup<HalfedgeDS,Polyhedron,Kernel,InputIterator> soup(begin,end);
		output.delegate(soup);
	}
};

#endif // _MAKE_QUAD_SOUP_
