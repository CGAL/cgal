#ifndef CGAL_SDG_LINF_LBIS_2_H
#define CGAL_SDG_LINF_LBIS_2_H

#include <CGAL/Segment_Delaunay_graph_Linf_2.h>

/*Sandeep for debugging
#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif*/
/*
  Conventions:
  ------------
  1. we treat segments as open; the endpoints are separate objects
  2. a segment of length zero is treated as a point
  3. a point is deleted only if it has no segment adjacent to it
  4. when a segment is deleted it's endpoints are not deleted
  5. the user can force the deletion of endpoints; this is only done
     if condition 3 is met.
  6. when objects are written to a stream we distinguish between
     points and segments; points start by a 'p' and segments by an 's'.
*/


namespace CGAL {

template<class Gt,
         class ST = Segment_Delaunay_graph_storage_traits_2<Gt>,
         class D_S = Triangulation_data_structure_2 < 
                Segment_Delaunay_graph_vertex_base_2<ST>,
                Segment_Delaunay_graph_face_base_2<Gt> >,
         class LTag = Tag_false >
class SDG_Linf_lbis_2
  : private Segment_Delaunay_graph_Linf_2<
          Gt, ST, D_S, LTag > 
{
};


template<class Gt, class D_S, class LTag>
std::istream& operator>>(std::istream& is,
			 SDG_Linf_lbis_2<Gt,D_S,LTag>& sdg)
{
  sdg.file_input(is);
  return is;
}

template<class Gt, class D_S, class LTag>
std::ostream& operator<<(std::ostream& os,
			 const SDG_Linf_lbis_2<Gt,D_S,LTag>& sdg)
{
  sdg.file_output(os);
  return os;
}

} //namespace CGAL


#endif // CGAL_SDG_LINF_LBIS_2_H
