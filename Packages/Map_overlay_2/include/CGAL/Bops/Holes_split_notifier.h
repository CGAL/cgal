#ifndef CGAL_HOLES_SPLIT_NOTIFIER_H
#define CGAL_HOLES_SPLIT_NOTIFIER_H   


CGAL_BEGIN_NAMESPACE         

template <class Planar_map_>
class Holes_split_notifier : public Planar_map_::Change_notification
{
public:
  typedef Planar_map_                                     Planar_map;
  typedef typename Planar_map::Halfedge_handle            Halfedge_handle;
  typedef typename Planar_map::Halfedge_const_handle      Halfedge_const_handle;
  
  typedef typename Planar_map::Traits                    Traits;
  typedef typename Traits::Point                          Point;
  typedef typename Traits::X_curve                        X_curve;
  
  Holes_split_notifier() : decomposing_(false) {}
  
  virtual ~Holes_split_notifier() {}
  
  // We assume the curves comes from the original triangles are oriented counter clockwise.
  void add_edge(const X_curve& cv, 
                Halfedge_handle e, 
                bool left_to_right, 
                bool overlap = false)
  { 
    if (decomposing_){
      e->set_decomposing(true);
      e->twin()->set_decomposing(true);
    }
  }
 
  
  /***************************************** new functions **********************************************************/
  
  
  void  set_decomposing_edge(bool decomposing) 
  { 
    decomposing_ = decomposing; 
  }

  bool  is_decomposing_edge() const { 
    return decomposing_; 
  }
  
private:
  bool  decomposing_;
};

CGAL_END_NAMESPACE

#endif




