#include <CGAL/Cartesian.h>
#include <CGAL/basic.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Viewer_stream.h>
#include <CGAL/facet_object.h>

typedef CGAL::Cartesian<double> rep_t;
typedef CGAL::Triangulation_geom_traits_3<rep_t>  traits_3;
typedef CGAL::Triangulation_vertex_base_3<traits_3>     Vb ;
typedef CGAL::Triangulation_cell_base_3<traits_3>       Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> TDS3 ;
typedef CGAL::Delaunay_triangulation_3<traits_3,TDS3> Delaunay_3;


CGAL::list_vertices get_vertices(const Delaunay_3 &fo)
{
  Delaunay_3::Vertex_iterator vit;
  CGAL::list_vertices res;
  CGAL::vertex v(3);
  for (vit=fo.finite_vertices_begin(); vit!=fo.vertices_end();vit++) {
    v[0]=CGAL::to_double(vit->point().x());
    v[1]=CGAL::to_double(vit->point().y());
    v[2]=CGAL::to_double(vit->point().z());
    res.push_back(v);
  }
return res;
}


CGAL::list_edges get_edges(const Delaunay_3 &fo)
{
  Delaunay_3::Edge_iterator it;
  CGAL::list_edges res;
  CGAL::edge ed(2);
  CGAL::vertex vx(3);
  Delaunay_3::Cell_handle f;
  int n1, n2;
  Delaunay_3::Vertex_handle v1, v2;
  for (it=fo.finite_edges_begin(); it!=fo.edges_end();it++) {
	f = (*it).first;
	n1 = (*it).second;
       	n2 = (*it).third;
	v1 = f->vertex(n1);
	v2 = f->vertex(n2);
    vx[0]=CGAL::to_double(v1->point().x());
    vx[1]=CGAL::to_double(v1->point().y());
    vx[2]=CGAL::to_double(v1->point().z());
    ed[0]=vx;
    vx[0]=CGAL::to_double(v2->point().x());
    vx[1]=CGAL::to_double(v2->point().y());
    vx[2]=CGAL::to_double(v2->point().z());
    ed[1]=vx;
    res.push_back(ed);
  }
return res;
}

CGAL::list_facets get_facets(const Delaunay_3 &fo)
{
  Delaunay_3::Facet_iterator it;
  CGAL::list_facets res;
  CGAL::facet fa;
  CGAL::vertex v(3);
  for (it=fo.finite_facets_begin(); it!=fo.facets_end();it++) {
    fa.clear();
    v[0]=CGAL::to_double((((*it).first)->vertex(((*it).second +1)
						&3))->point().x());

    v[1]=CGAL::to_double((((*it).first)->vertex(((*it).second +1)
						&3))->point().y());

    v[2]=CGAL::to_double((((*it).first)->vertex(((*it).second +1)
						&3))->point().z());

    fa.push_back(v);
    v[0]=CGAL::to_double((((*it).first)->vertex(((*it).second +2)
						&3))->point().x());

    v[1]=CGAL::to_double((((*it).first)->vertex(((*it).second +2)
						&3))->point().y());

    v[2]=CGAL::to_double((((*it).first)->vertex(((*it).second +2)
						&3))->point().z());

    fa.push_back(v);
    v[0]=CGAL::to_double((((*it).first)->vertex(((*it).second +3)
						&3))->point().x());

    v[1]=CGAL::to_double((((*it).first)->vertex(((*it).second +3)
						&3))->point().y());

    v[2]=CGAL::to_double((((*it).first)->vertex(((*it).second +3)
						&3))->point().z());

    fa.push_back(v);
    res.push_back(fa);
  }
return res;
}
