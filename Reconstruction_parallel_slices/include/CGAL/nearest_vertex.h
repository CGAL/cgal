  
namespace CGAL {


template < typename Delaunay>
typename Delaunay::Vertex_handle
  nearest_vertex(const Delaunay& dt, const typename Delaunay::Point& p, typename Delaunay::Face_handle f)
{
  CGAL_triangulation_precondition(dt.dimension() == 2);
  f = dt.locate(p,f);

  typename Delaunay::Geom_traits::Compare_distance_2 
    compare_distance =  dt.geom_traits().compare_distance_2_object();
  typename Delaunay::Vertex_handle nn =  !dt.is_infinite(f->vertex(0)) ? f->vertex(0):f->vertex(1);
  if ( !dt.is_infinite(f->vertex(1)) && compare_distance(p,
					    f->vertex(1)->point(),
					    nn->point()) == SMALLER) 
    nn=f->vertex(1);
  if ( !dt.is_infinite(f->vertex(2)) && compare_distance(p,
					    f->vertex(2)->point(), 
					    nn->point()) == SMALLER) 
    nn=f->vertex(2);
       
  look_nearest_neighbor(dt, p,f,0,nn);
  look_nearest_neighbor(dt, p,f,1,nn);
  look_nearest_neighbor(dt, p,f,2,nn);
  return nn;
}


template < typename Delaunay>
void
look_nearest_neighbor(const Delaunay& dt,
                      const typename Delaunay::Point& p,
                      typename Delaunay::Face_handle f,
		      int i,
		       typename Delaunay::Vertex_handle& nn)
{
   typename Delaunay::Face_handle  ni=f->neighbor(i);
  if ( ON_POSITIVE_SIDE != dt.side_of_oriented_circle(ni,p,true) ) return;

  typename Delaunay::Geom_traits::Compare_distance_2 
    compare_distance =  dt.geom_traits().compare_distance_2_object();
  i = ni->index(f);
  if ( !dt.is_infinite(ni->vertex(i)) &&
       compare_distance(p, 
	      ni->vertex(i)->point(),
	      nn->point())  == SMALLER)  nn=ni->vertex(i);
    
  // recursive exploration of triangles whose circumcircle contains p
  look_nearest_neighbor(dt, p, ni, Delaunay::ccw(i), nn);
  look_nearest_neighbor(dt, p, ni, Delaunay::cw(i), nn);
} 

}
