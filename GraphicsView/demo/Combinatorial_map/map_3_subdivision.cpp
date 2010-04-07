#include "typedefs.h"

#define PI 3.1415926535897932

// Smoth a vertex depending on the vertices of its incident face.
class Smooth_old_vertex
{
public:
  /**  Contructor.
   * @param amap is the map to smooth
   * @param amark is a mark designing old darts (i.e. darts not created during
   *        the triangulation step)
   */
  Smooth_old_vertex(Map &amap, unsigned int amark) : mmap(amap)
  {}
  
  Vertex operator()( Vertex& v) const
  {
    Dart_handle d = v.dart();
    CGAL_assertion(d!=NULL);

    int degree=0;
    bool open = false;

    Map::Edge_iterator_of_vertex it(mmap, d);
    for ( ; it!=mmap.edge_iterator_of_vertex_end(d); ++it )
    {
      ++degree;
      if ( it->is_free(2) ) open = true;
    }        

    if ( open ) return v;

    Map::FT alpha = (4.0f - 2.0f *
                     (Map::FT)cos( 2.0f * PI / (Map::FT)degree)) / 9.0f;
    Map::Vector vec = (v.point() - CGAL::ORIGIN) * ( 1.0f - alpha);

    for (it.rewind(); it!=mmap.edge_iterator_of_vertex_end(d); ++it)
    {
      CGAL_assertion(!it->is_free(2));

        vec = vec + (it->second_vertex()->point() - CGAL::ORIGIN)
        * alpha / degree;
    }
    
    Vertex res(CGAL::ORIGIN + vec);
    res.set_dart(d);
    return res;
  }
private:
  Map& mmap;
};

// Flip an edge, work in 2D and in 3D.
Dart_handle flip_edge(Map &m, Dart_handle d)
{
  CGAL_assertion( d!=NULL && !d->is_free(2) );

  if ( !m.can_remove(d,1) ) return NULL;

  Dart_handle d2 = d->beta(1,1);
  remove_edge_3(m, d);
  insert_edge_3(m, d2, d2->beta(1, 1));

  return d2->beta(0);
}

// Subdivide each face of the map by using sqrt(3)-subdivision.
void subdivide_map_3(Map& m)
{
  if (m.size_of_darts() == 0)
    return;

  unsigned int mark    = m.get_new_mark();
  unsigned int treated = m.get_new_mark();
  m.negate_mask_mark(mark); // All the old darts are marked in O(1).

  // 1) We smoth the old vertices.
  std::vector<Vertex> vertices; // smooth the old vertices
  vertices.reserve(m.size_of_vertices()); // get intermediate space
  std::transform( m.vertices_begin(), m.vertices_end(),
                  std::back_inserter(vertices), Smooth_old_vertex(m,mark));

  // 2) We subdivide each face.
  m.negate_mask_mark(treated); // All the darts are marked in O(1).
  unsigned int nb=0;
  for ( Map::All_darts_iterator it(m.darts_begin());
        m.number_of_marked_darts(treated)>0; ++it )
  {
    ++nb;
    if ( m.is_marked(it, treated) )
    {
      // We unmark the darts of the face to process only once dart/face.
      CGAL::unmark_orbit<Map>(m, it, Map::FACE_ORBIT, treated);
      // We triangulate the face.
      CGAL::triangulate_face_3(m, it);
    }
  }

  CGAL_assertion( m.is_whole_map_unmarked(treated) );
  CGAL_assertion( m.is_valid() );
  m.free_mark(treated);

  // 3) We update the coordinates of old vertices.
  for(std::vector<Vertex>::iterator vit=vertices.begin();
  vit!=vertices.end();++vit)
  {
    CGAL_assertion(vit->dart()!=NULL);
    CGAL_assertion(vit->dart()->vertex()!=NULL);
    vit->dart()->vertex()->set_point(vit->point());
  }

  // 4) We flip all the old edges.
  m.negate_mask_mark(mark); // Now only new darts are marked.
  Dart_handle d2 = NULL;
  for (Map::All_darts_iterator it(m.darts_begin()); it != m.darts_end();)
  {
    d2 = it++;
    CGAL_assertion(d2!=NULL);
     if (!m.is_marked(d2, mark))   // This is an old dart.
    {
      // We process only the last dart of a same edge.
      if (!d2->is_free(2) && (d2->beta(2,3)==d2->beta(3,2)))
      {
        if ( m.is_marked(d2->beta(2), mark) &&
             (d2->is_free(3) ||
              (m.is_marked(d2->beta(3), mark) &&
                m.is_marked(d2->beta(2,3), mark))) )
        {
          m.negate_mask_mark(mark); // thus new darts will be marked
          flip_edge(m, d2);
          m.negate_mask_mark(mark);
        }
        else
          m.set_mark(d2, mark);
      }
      else
        m.set_mark(d2, mark);
    }
  }
  CGAL_assertion( m.is_whole_map_marked(mark) );
  m.free_mark(mark);

  CGAL_postcondition(m.is_valid());
}
