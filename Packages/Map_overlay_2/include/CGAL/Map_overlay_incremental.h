#ifndef CGAL_MAP_OVERLAY_INCREMENTAL_H
#define CGAL_MAP_OVERLAY_INCREMENTAL_H

#include <vector>
#include <list>

#ifndef CGAL_MAP_OVERLAY_BASE_H
#include <CGAL/Map_overlay_base.h>
#endif

#ifndef CGAL_MAP_OVERLAY_DEFAULT_NOTIFIER_H
#include <CGAL/Map_overlay_default_notifier.h>
#endif

//#include <CGAL/IO/leda_window.h>  //used for visualization -
//should be placed after the other CGAL files.

CGAL_BEGIN_NAMESPACE

template  <class Arrangement_, 
  class Map_overlay_ChangeNotification_ = Map_overlay_default_notifier<Arrangement_> > 
class  Map_overlay_incremental : 
  public Map_overlay_base<Arrangement_, Map_overlay_ChangeNotification_> {
public:  
  
  typedef Arrangement_                                          Arrangement;
  typedef Map_overlay_ChangeNotification_                       Map_overlay_ChangeNotification;
  //typedef typename Arrangement::Planar_map                      PM;
  //typedef typename Arrangement::Pmwx                          Pmwx;
  typedef typename Arrangement::Vertex                          Vertex;
  typedef typename Arrangement::Face                            Face;
  typedef typename Arrangement::Halfedge                        Halfedge;
  typedef typename Arrangement::Vertex_handle                   Vertex_handle;
  typedef typename Arrangement::Halfedge_handle                 Halfedge_handle;
  typedef typename Arrangement::Face_handle                     Face_handle;
  typedef typename Arrangement::Vertex_const_handle             Vertex_const_handle;
  typedef typename Arrangement::Halfedge_const_handle           Halfedge_const_handle;
  typedef typename Arrangement::Face_const_handle               Face_const_handle;
  typedef typename Arrangement::Vertex_iterator                 Vertex_iterator;
  typedef typename Arrangement::Halfedge_iterator               Halfedge_iterator;
  typedef typename Arrangement::Face_iterator                   Face_iterator;
  typedef typename Arrangement::Vertex_const_iterator           Vertex_const_iterator;
  typedef typename Arrangement::Halfedge_const_iterator         Halfedge_const_iterator;
  typedef typename Arrangement::Face_const_iterator             Face_const_iterator;
  typedef typename Arrangement::Ccb_halfedge_circulator         Ccb_halfedge_circulator;
  typedef typename Arrangement::Ccb_halfedge_const_circulator   Ccb_halfedge_const_circulator;
  typedef typename Arrangement::Holes_iterator                  Holes_iterator;
  typedef typename Arrangement::Holes_const_iterator            Holes_const_iterator;
  typedef typename Arrangement::Locate_type                     Locate_type;

  typedef typename  Arrangement::Traits                    Traits;
  typedef typename  Traits::X_curve                        X_curve;
  typedef typename  Traits::Point                          Point;

  /*typedef typename PM::Vertex_handle                            Pm_vertex_handle;
    typedef typename PM::Halfedge_handle                          Pm_halfedge_handle;
    typedef typename PM::Face_handle                              Pm_face_handle;
    typedef typename PM::Vertex_const_handle                      Pm_vertex_const_handle;
    typedef typename PM::Halfedge_const_handle                    Pm_halfedge_const_handle;
    typedef typename PM::Face_const_handle                        Pm_face_const_handle;*/
  
private:
  typedef Map_overlay_incremental<Arrangement, Map_overlay_ChangeNotification>       Self;
public:

  Map_overlay_incremental() {}

  ~Map_overlay_incremental() {}
  
  void  map_overlay(const Arrangement &a1, 
                    const Arrangement &a2, 
                    Map_overlay_ChangeNotification *change_notf, 
                    Arrangement &result)
  {
    Halfedge_const_iterator  halfedge_iter;
    
    // Creaing the new arrangement induced by both subdivisions. 
    // inserting the curves of every halfedge one by one.
    
    for (halfedge_iter = a1.halfedges_begin(); 
         halfedge_iter != a1.halfedges_end(); 
         ++halfedge_iter, ++halfedge_iter){
      change_notf->set_curve_attributes(halfedge_iter->curve(), 
                                        halfedge_iter, 
                                        true);
      result.insert(halfedge_iter->curve(), change_notf);
    }
    
    for (halfedge_iter = a2.halfedges_begin(); 
         halfedge_iter != a2.halfedges_end(); 
         ++halfedge_iter, ++halfedge_iter){
      change_notf->set_curve_attributes(halfedge_iter->curve(), 
                                        halfedge_iter, 
                                        false); 
      result.insert(halfedge_iter->curve(), change_notf);
    }
    
    //cout<<"After creating Arr\n";
    change_notf->update_all_faces(result);
  }

};

CGAL_END_NAMESPACE

#endif











