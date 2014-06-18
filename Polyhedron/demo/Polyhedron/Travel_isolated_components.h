#ifndef TRAVEL_ISOLATED_COMPONENTS_H
#define TRAVEL_ISOLATED_COMPONENTS_H

#include <boost/optional.hpp>
#include <vector>
#include "One_ring_iterators.h"
#include <CGAL/Default.h>
class Travel_isolated_components {
public:
  // for transform iterator
  template<class HandleType>
  struct Get_handle {
    typedef HandleType result_type;
    template<class Iterator>
    result_type operator()(Iterator it) const
    { return it; }
  };

  // to be used in get_minimum_isolated_component function
  struct Minimum_visitor
  {
    template<class HandleType>
    void operator()(const std::vector<HandleType>& C) {
      if(!minimum) { minimum = C.size(); }
      else         { minimum = (std::min)(*minimum, C.size()); }
    }

    boost::optional<std::size_t> minimum;
  };

  // to be used in select_isolated_components function
  template<class OutputIterator>
  struct Selection_visitor
  {
    Selection_visitor(std::size_t threshold_size, OutputIterator out) 
      : threshold_size(threshold_size), out(out), any_inserted(false) { }

    template<class HandleType>
    void operator()(const std::vector<HandleType>& C) {
      if(C.size() <= threshold_size) {
        any_inserted = true;
        out = std::copy(C.begin(), C.end(), out);
      }
      else {
        minimum_visitor(C);
      }
    }

    std::size_t     threshold_size;
    OutputIterator  out;
    bool            any_inserted;
    Minimum_visitor minimum_visitor; // hold minimum of NOT inserted components
  };

  // NOTE: prior to call this function, id fields should be updated
  template<class HandleType, class InputIterator, class IsSelected, class Visitor>
  void travel(InputIterator begin, 
              InputIterator end, 
              std::size_t size, 
              const IsSelected& selection, 
              Visitor& visitor)
  {
    std::vector<bool> mark(size, false);

    for(; begin != end; ++begin) 
    {
      HandleType h = begin;

      if(mark[h->id()] || selection.is_selected(h)) { continue; }

      std::vector<HandleType> C;
      C.push_back(h);
      mark[h->id()] = true;
      std::size_t current_index = 0;

      bool neigh_to_selection = false;
      while(current_index < C.size()) {
        HandleType current = C[current_index++];

        for(One_ring_iterator<HandleType> circ(current); circ; ++circ)
        {
          HandleType nv = circ;
          neigh_to_selection |= selection.is_selected(nv);
          if(!mark[nv->id()] && !selection.is_selected(nv)) {
            mark[nv->id()] = true; 
            C.push_back(nv);
          }
        }
      }
      if(neigh_to_selection) { visitor(C); }
    }
  }
};
#endif
