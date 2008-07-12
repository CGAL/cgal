#ifndef RENDER
#define RENDER

#include <GL/glu.h>

// forward declaration
template <class Facet, class Kernel>
typename Kernel::Vector_3 compute_normal(Facet& f);

template <class Polyhedron>
void gl_render(Polyhedron& polyhedron)
{
	typedef typename Polyhedron::Traits Kernel;
	typedef typename Kernel::Point_3 Point;
	typedef typename Kernel::Vector_3 Vector;

	typedef typename Polyhedron::Halfedge_handle                    Halfedge_handle;
	typedef typename Polyhedron::Facet_handle                       Facet_handle;
	typedef typename Polyhedron::Vertex_handle                      Vertex_handle;

	typedef typename Polyhedron::Vertex_iterator                    Vertex_iterator;
	typedef typename Polyhedron::Facet_iterator                     Facet_iterator;
	typedef typename Polyhedron::Halfedge_iterator                  Halfedge_iterator;
	typedef typename Polyhedron::Edge_iterator                      Edge_iterator;

	typedef typename Polyhedron::Halfedge_around_vertex_circulator  HV_circulator;
	typedef typename Polyhedron::Halfedge_around_facet_circulator   HF_circulator;

	Facet_iterator f;
	for(f = polyhedron.facets_begin();
	    f != polyhedron.facets_end();
	    f++)
	{
		::glBegin(GL_POLYGON);

		// compute normal
		Vector n = compute_normal<typename Polyhedron::Facet,Kernel>(*f);
		::glNormal3d(n.x(),n.y(),n.z());

		// revolve around current face to get vertices
		HF_circulator he = f->facet_begin();
		do
		{
			const Point& p = he->vertex()->point();
			::glVertex3d(p.x(),p.y(),p.z());
		}
		while(++he != f->facet_begin());
		::glEnd();
	}
}

template <class Facet, class Kernel>
typename Kernel::Vector_3 compute_normal(Facet& f)
{
	typedef typename Kernel::Point_3 Point;
	typedef typename Kernel::Vector_3 Vector;
	typedef typename Facet::Halfedge_around_facet_circulator HF_circulator;
	Vector normal = CGAL::NULL_VECTOR;
	HF_circulator he = f.facet_begin();
	do
	{
		const Point& prev = he->prev()->vertex()->point();
		const Point& curr = he->vertex()->point();
		const Point& next = he->next()->vertex()->point();
		Vector n = CGAL::cross_product(next-curr,prev-curr);
		normal = normal + (n / std::sqrt(n*n));
	}
	while (++he != f.facet_begin());
	return normal / std::sqrt(normal * normal);
}

#endif // RENDER




