#ifndef CGAL_MAP_OVERLAY_NAIVE_H
#define CGAL_MAP_OVERLAY_NAIVE_H

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
class  Map_overlay_naive : 
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
  typedef Map_overlay_naive<Arrangement, Map_overlay_ChangeNotification>       Self;
public:

  Map_overlay_naive() {}

  void  map_overlay(const Arrangement &a1, 
                    const Arrangement &a2, 
                    Map_overlay_ChangeNotification *pmwx_change_notf, 
                    Arrangement &result)
  {
    /*typedef typename Arrangement::Vertex_handle             Vertex_handle;
      typedef typename Arrangement::Halfedge_handle           Halfedge_handle;
      typedef typename Arrangement::Face_handle               Face_handle;
      typedef typename Arrangement::Vertex_const_handle       Vertex_const_handle;
      typedef typename Arrangement::Halfedge_const_handle     Halfedge_const_handle;
      typedef typename Arrangement::Face_const_handle         Face_const_handle;
    
      typedef typename  Arrangement::Vertex                    Vertex;
      typedef typename  Arrangement::Halfedge                  Halfedge;
      typedef typename  Arrangement::Face                      Face;
      
      typedef typename  Arrangement::Halfedge_const_iterator   Halfedge_const_iterator; 
      typedef typename  Arrangement::Face_iterator             Face_iterator;
      typedef typename  Arrangement::Face_const_iterator       Face_const_iterator;
      typedef typename  Arrangement::Ccb_halfedge_circulator   Ccb_halfedge_circulator;
      typedef typename  Arrangement::Planar_map                PM;
      typedef typename  Arrangement::Locate_type               Locate_type;
      typedef typename  Arrangement::Traits                    Traits;
      typedef typename  Traits::X_curve                        X_curve;
      typedef typename  Traits::Point                          Point;*/
    
    
    Halfedge_const_iterator  half_edge_iter;
    
    // debuging only!
    //if (pmwx_change_notf == NULL)
    //  std::cout<<"Notifier is NULL" << std::endl;

    // Creaing the new arrangement of both subdivisions. inserting the curevs of the halfedges (already calculated fomr the given subdivisions 
    for (half_edge_iter = a1.halfedges_begin(); 
         half_edge_iter != a1.halfedges_end(); 
         ++half_edge_iter, ++half_edge_iter){
      pmwx_change_notf->set_curve_attributes(half_edge_iter->curve(), 
                                             half_edge_iter, 
                                             true);
      result.insert(half_edge_iter->curve(), pmwx_change_notf);
      
      //result.insert(half_edge_iter->curve(), NULL); //debug only!
    }
    
    for (half_edge_iter = a2.halfedges_begin(); 
         half_edge_iter != a2.halfedges_end(); 
         ++half_edge_iter, ++half_edge_iter){
      pmwx_change_notf->set_curve_attributes(half_edge_iter->curve(), 
                                             half_edge_iter, 
                                             false); 
      result.insert(half_edge_iter->curve(), pmwx_change_notf);
      
      //result.insert(half_edge_iter->curve(), NULL);  //debug only!
    }
    
    //cout<<"After creating Arr\n";
  }

};

CGAL_END_NAMESPACE

#endif











