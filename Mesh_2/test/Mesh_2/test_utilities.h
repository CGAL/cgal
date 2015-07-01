
#ifndef CGAL_TEST_UTILITIES_2_H
#define CGAL_TEST_UTILITIES_2_H


template <class CTr>
typename CTr::size_type number_of_constrained_edges(const CTr& tr)
{
  typename CTr::size_type nedges = 0;
  for(typename CTr::Finite_edges_iterator eit = tr.finite_edges_begin();
      eit != tr.finite_edges_end();
      ++eit)
    if(tr.is_constrained(*eit))
      ++nedges;
  return nedges;
}


#endif //CGAL_TEST_UTILITIES_2_H
