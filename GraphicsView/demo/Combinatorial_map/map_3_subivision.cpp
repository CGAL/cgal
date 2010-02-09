#include "typedefs.h"

#define PI 3.1415926535897932

// Smoth a vertex depending on the vertices of its incident face.
class Smooth_old_vertex
{
public:
  Smooth_old_vertex(Map &amap) : m(amap)
  {}
  
  Point_3 operator()( Vertex& v) const
  {
    Dart_handle d = v.dart();
    CGAL_assertion(d!=NULL);
    
    int degree=0;
    bool open = false;
    
    Map::Edge_iterator_of_vertex it(m, d);   
    for ( ; it!=m.edge_iterator_of_vertex_end(d); ++it )
      {
	// std::cout<<"Dart : "<<(&**it)<<std::endl;
	++degree;
	if (it->is_free(2) ) open = true;
      }
    
    if (open || degree & 1 )
      {
	//	std::cout<<"Smooth_old_vertex : no smooth "<<v.point()<<std::endl;
	return v.point();
      }

    degree /= 2;
    
    Map::FT alpha = (4.0f - 2.0f *
		     (Map::FT)cos( 2.0f * PI / (Map::FT)degree)) / 9.0f;
    Map::Vector vec = (v.point() - CGAL::ORIGIN) * ( 1.0f - alpha);


    //    std::cout<<"Degre:"<<degree<<std::endl<<"Point : "<<v.point()<<std::endl;	

    for (it.rewind(); it!=m.edge_iterator_of_vertex_end(d); ++it,++it)
      {
	CGAL_assertion(!it->is_free(2));
	
	//	std::cout<<"  pt : "<<it->second_vertex()->point()<<std::endl;
	
	vec = vec + (it->second_vertex()->point() - CGAL::ORIGIN)
	  * alpha / degree;
      }
    
    //    std::cout<<"Smooth_old_vertex : "<<CGAL::ORIGIN + vec<<std::endl;
    
    return (CGAL::ORIGIN + vec);
  }
private:
  Map& m;
};

// Flip an edge, work in 2D and in 3D.
Dart_handle flip_edge(Map &m, Dart_handle d)
{
   CGAL_assertion( d!=NULL && !d->is_free(2) );
      
   if ( !can_remove(m,d,1) ) return NULL;
   
   Dart_handle d2 = d->beta(1,1);
   remove_edge_3(m, d);
   insert_edge_3(m, d2, d2->beta(1, 1));
   return d2->beta(0);
}

// Subdivide each face of the map by using sqrt(3)-subdivision.
void subdivide_map_3(Map& m)
{
   if (m.number_of_darts() == 0)
      return;
   
   m.display_characteristics(std::cout)<<std::endl<<"Nb vertices:"<< m.size_of_vertices()<<std::endl;
   
     // We use that new vertices/halfedges/facets are appended at the end.
   std::size_t nv = m.size_of_vertices();
   Map::Vertex_iterator last_v = m.vertices_end();
   -- last_v;  // the last of the old vertices

   unsigned int mark    = m.get_new_mark();
   unsigned int treated = m.get_new_mark();
   m.negate_mask_mark(mark); // All the old darts are marked in O(1).

   // 1) We subdivide each face.
   for ( Map::All_darts_iterator it(m.darts_begin()); 
         it!=m.darts_end(); ++it )
     {
       if ( m.is_marked(it, mark) && !m.is_marked(it, treated) )
	 {
	   // We mark the darts of the face to process only once dart/face. 
	   CGAL::mark_orbit<Map>(m, it, Map::FACE_ORBIT, treated);
	   // We triangulate the face.
	   CGAL::triangulate_face_3(m, it);
	 }
       else m.set_mark(it, treated);   
     }
   
   CGAL_assertion( m.is_whole_map_marked(treated) );
   m.free_mark(treated);
   CGAL_assertion(m.is_valid());

   // 2) We smoth the old vertices.
   std::vector<Point_3> pts; // smooth the old vertices
   pts.reserve(nv);          // get intermediate space for the new points
   ++ last_v;                // make it the past-the-end position again
   std::transform( m.vertices_begin(), last_v, std::back_inserter(pts),
		   Smooth_old_vertex(m));
   std::copy( pts.begin(), pts.end(), m.vertices_begin());

   /*   std::vector<Point_3>::iterator pts_it=pts.begin();   
   for(Map::Vertex_iterator it=m.vertices_begin();it!=last_v;++it)
     {
       it->set_point(*pts_it);
       ++pts_it;
       }*/   
   
   // 3) We flip all the old edges.   
   m.negate_mask_mark(mark); // Now only new darts are marked.
   Dart_handle d2 = NULL;
   for (Map::All_darts_iterator it(m.darts_begin()); it != m.darts_end();)
   {
      d2 = it++;
      if (!m.is_marked(d2, mark))   // This is an old dart.
      {
         // We process only the last dart of a same edge.
	if (!d2->is_free(2) && m.is_marked(d2->beta(2), mark) &&
	    (d2->is_free(3) ||
	     (m.is_marked(d2->beta(3), mark) && m.is_marked(d2->beta(2,3), mark))))
	  {
            m.negate_mask_mark(mark); // thus new darts will be marked
            flip_edge(m, d2);
            m.negate_mask_mark(mark);
	  }
	else
	  m.set_mark(d2, mark);
      }
   }
   CGAL_assertion( m.is_whole_map_marked(mark) );
   m.free_mark(mark);

   CGAL_postcondition(m.is_valid());
}
