#ifndef CGAL_AFSR_WRITE_TRIPLE_INDICES_H
#define CGAL_AFSR_WRITE_TRIPLE_INDICES_H

#include <CGAL/tuple.h>

namespace CGAL {

template <class Kernel, class Triangulation, class Filter>
class Advancing_front_surface_reconstruction;



template <class OutputIterator, class Kernel, class Triangulation, class Filter>
OutputIterator
  write_triple_indices(OutputIterator out, const Advancing_front_surface_reconstruction<Kernel,Triangulation,Filter>& S)
{ 
  typedef Advancing_front_surface_reconstruction<Kernel,Triangulation,Filter> Surface;
  typedef typename Surface::TDS_2 TDS_2;
  typedef typename TDS_2::Face_iterator Face_iterator;

  if(S.triangulation_3().dimension() < 3){
    std::cerr << "not 3D\n";
    return out;
  }
  const TDS_2& tds = S.triangulation_data_structure_2();

  for(Face_iterator fit = tds.faces_begin(); fit != tds.faces_end(); ++fit){

    if(fit->is_on_surface()){
      *out++ = CGAL::cpp11::tuple<std::size_t,std::size_t,std::size_t>(fit->vertex(0)->vertex_3()->id(),
                                                                       fit->vertex(1)->vertex_3()->id(),
                                                                       fit->vertex(2)->vertex_3()->id());
    }
  }
    return out;
}

} 


#endif
