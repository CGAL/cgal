#ifndef CGAL_APOLLONIUS_GRAPH_EUCLIDEAN_TRAITS_WRAPPER_2_H
#define CGAL_APOLLONIUS_GRAPH_EUCLIDEAN_TRAITS_WRAPPER_2_H

CGAL_BEGIN_NAMESPACE


template<class Gt_base>
class Apollonius_graph_gt_wrapper : public Gt_base
{
public:
  typedef typename Gt_base::Rep::Point_2    Triangle_2;

  Apollonius_graph_gt_wrapper() {}
  Apollonius_graph_gt_wrapper(const Apollonius_graph_gt_wrapper&) {}
  Apollonius_graph_gt_wrapper
  operator=(const Apollonius_graph_gt_wrapper&) {
    return (*this);
  }

  Apollonius_graph_gt_wrapper(const Gt_base&) {}
};




CGAL_END_NAMESPACE


#endif // CGAL_APOLLONIUS_GRAPH_EUCLIDEAN_TRAITS_WRAPPER_2_H
