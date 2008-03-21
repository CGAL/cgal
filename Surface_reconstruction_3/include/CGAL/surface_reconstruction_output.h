#ifndef CGAL_SURFACE_RECONSTRUCTION_OUTPUT_H
#define CGAL_SURFACE_RECONSTRUCTION_OUTPUT_H

CGAL_BEGIN_NAMESPACE


template <class C2t3, class Triangle>
void
output_surface_facets(std::list<Triangle>& triangles,
					            const C2t3& c2t3)
{
    typedef typename C2t3::Triangulation Tr;
    typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Point Point;

    const Tr& tr = c2t3.triangulation();

    for(Finite_facets_iterator fit = tr.finite_facets_begin();
        fit != tr.finite_facets_end();
		    ++fit)
	{
        if((*fit).first->is_facet_on_surface((*fit).second) == true)
        {
            std::vector<Point> points(3);
            unsigned int index = 0;
            for(int i=0; i<4; i++)
                if(i != (*fit).second)
		            points[index++] = (*fit).first->vertex(i)->point();
            triangles.push_back(Triangle(points[0],points[1],points[2]));
        }
	}
}


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_RECONSTRUCTION_OUTPUT_H
