// Author     : Nader Salman


// an object contais the constrained triangle soup

#ifndef _Gyroviz_constrained_triangle_soup_
#define _Gyroviz_constrained_triangle_soup_	

#include "Gyroviz_triangle_with_cam.h"

template < class Triangulation >
class Gyroviz_constrained_triangle_soup : public std::vector<Gyroviz_triangle_with_cam <Triangulation> >
{

private:
	// Vertex_handle Triangle_3
	typedef typename Triangulation::Vertex_handle      Vertex_handle;
	typedef typename Gyroviz_triangle_with_cam<Triangulation>   Gyroviz_triangle_with_cam;
	/*std::vector<Gyroviz_triangle_with_cam> vector_cam;*/

public:

	// constructor
	Gyroviz_constrained_triangle_soup(){}
	Gyroviz_constrained_triangle_soup(std::vector<Gyroviz_triangle_with_cam> vc):this(vc){}

	// accessors
	//const	std::vector<Gyroviz_triangle_with_cam>&	get_vector_cam()    const { return	vector_cam; }

	// 3D projection of the tracked 2D constrained triangulation
	void gl_draw_soup_triangles(const unsigned char r, const unsigned char g,
		const unsigned char b){

			::glColor3ub(r,g,b);
			::glBegin(GL_TRIANGLES);

			int test_size = this->size();

			for(int i=0; i<this->size(); ++i)
			{
				Point_3 p1 = (*this)[i].get_triangle3()[0];
				Point_3 p2 = (*this)[i].get_triangle3()[1];
				Point_3 p3 = (*this)[i].get_triangle3()[2];
				::glVertex3d(p1.x(), p1.y(), p1.z());
				::glVertex3d(p2.x(), p2.y(), p2.z());
				::glVertex3d(p3.x(), p3.y(), p3.z());
			} 
			::glEnd();
	} 
};

#endif // _Gyroviz_constrained_triangle_soup_