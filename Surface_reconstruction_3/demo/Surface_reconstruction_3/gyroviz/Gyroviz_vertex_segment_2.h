// Author     : Nader Salman


// an object used to store a segment_2 as a couple of vertices Source & Target
// will be used to insert constraints in the 2D Delaunay triangulation
//  function:   ct.insert_constraint ( Vertex_handle va, Vertex_handle vb)
//		           Inserts the line segment s whose endpoints are the vertices 
//                 va and vb as a constrained edge e. The triangles intersected
//			       by s are removed and new ones are created.

#ifndef _Gyroviz_vertex_segment_2_
#define _Gyroviz_vertex_segment_2_

#include <CGAL/Constrained_Delaunay_triangulation_2.h>

template < class Triangulation >
class	Gyroviz_vertex_segment_2
{

protected:

	typedef typename Triangulation::Vertex_handle      Vertex_handle;

	// Vertex_handles
	Vertex_handle source;
	Vertex_handle target;


public:
	Gyroviz_vertex_segment_2(){}
	Gyroviz_vertex_segment_2(Vertex_handle va, Vertex_handle vb):source(va),target(vb){}

	// accessors
	const	Vertex_handle	get_source() const { return	source; }
	const	Vertex_handle	get_target() const { return	target; }

	// modificators
	void	set_source(Vertex_handle va) { source=va; }
	void	set_target(Vertex_handle vb) { target=vb; }
};

#endif // _Gyroviz_vertex_segment_2_