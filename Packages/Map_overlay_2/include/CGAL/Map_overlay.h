#ifndef CGAL_MAP_OVERLAY_H
#define CGAL_MAP_OVERLAY_H

#include <vector>
#include <list>
#include <time.h>

#ifndef CGAL_MAP_OVERLAY_BASE_H
#include <CGAL/Map_overlay_base.h>
#endif

#ifndef CGAL_MAP_OVERLAY_SWEEP_H
#include <CGAL/Map_overlay_sweep.h>
#endif

#ifndef CGAL_MAP_OVERLAY_DEFAULT_NOTIFIER_H
#include <CGAL/Map_overlay_default_notifier.h>
#endif


CGAL_BEGIN_NAMESPACE

template  <class Arrangement_, 
           class Change_notification_ = Map_overlay_default_notifier<Arrangement_> > 
class  Map_overlay {
public:  
  
  typedef Arrangement_                                          Arrangement;
  typedef Change_notification_                                  Change_notification;
  typedef typename Arrangement::Planar_map                      Planar_map;
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
  //typedef typename Arrangement::Point_location_base             Point_location_base;

  typedef Map_overlay_base<Arrangement, Change_notification>   Map_ovl_base;
  typedef Map_overlay_sweep<Arrangement, Change_notification>  Map_ovl_sweep;
  
  typedef Pm_point_location_base<Planar_map> Point_location_base;
  
private:
  typedef Map_overlay<Arrangement, Change_notification>       Self;
  
public:

  Map_overlay() : sub_division1(0), sub_division2(0), 
    ovl_change_notf(new Change_notification),  
    ovl(new Map_ovl_sweep), 
    use_delete_notf(true), 
    use_delete_ovl(true) {
    
    /*ovl = new Map_ovl_sweep;
      use_delete_ovl = true;
      
      ovl_change_notf = new Map_overlay_change_notification;
      use_delete_notf = true;*/
  }

  Map_overlay (const Arrangement &arr) : arr_(arr), 
    sub_division1(0), sub_division2(0),  
    ovl_change_notf(new Change_notification),  
    ovl(new Map_ovl_sweep), 
    use_delete_notf(true), 
    use_delete_ovl(true)
  { 
    //copy_arr(a, arr);
    
    /*  ovl = new Map_ovl_sweep;
        use_delete_ovl = true;
    
        ovl_change_notf = new Map_overlay_change_notification;
        use_delete_notf = true;
        
        sub_division1 = sub_division2 = NULL;*/
  }
  
  Map_overlay (const Arrangement &arr, 
               Change_notification* pmwx_change_notf) : 
    arr_(arr), sub_division1(0), sub_division2(0), 
    ovl_change_notf(pmwx_change_notf), ovl(new Map_ovl_sweep), 
    use_delete_notf(false), use_delete_ovl(true)  
  { 
    // here we can't use copy Constructor since we have to update arr attributres due to the notifier.
    // An effeicient way doing this is sweeping the original subdivision while using the notifier.
    
    //Arrangement  empty_subdivision;
    //ovl->map_overlay(arr, empty_subdivision, ovl_change_notf, arr_);
  }
  
  Map_overlay (const Arrangement &arr, 
               Map_ovl_base *ovl_ptr) : 
    arr_(arr),  sub_division1(0), sub_division2(0), 
    ovl_change_notf(new Change_notification), ovl(ovl_ptr), 
    use_delete_notf(true), use_delete_ovl(false) {}

  Map_overlay (const Arrangement &arr,  
               Change_notification* pmwx_change_notf, 
               Map_ovl_base *ovl_ptr) : 
    arr_(arr), sub_division1(0), sub_division2(0), 
    ovl_change_notf(pmwx_change_notf), ovl(ovl_ptr),
    use_delete_notf(false), use_delete_ovl(false) {}

  Map_overlay (const Self &ovl1, const Self &ovl2) : 
    sub_division1(&ovl1), sub_division2(&ovl2),
    ovl_change_notf(new Change_notification(&(ovl1.subdivision()), 
                                            &(ovl2.subdivision()) )), 
    ovl(new Map_ovl_sweep), 
    use_delete_notf(true), use_delete_ovl(true)
  {
    //ovl = new Map_ovl_sweep;
    //use_delete_ovl = true;

    //ovl_change_notf = new Change_notification( &(ovl1.subdivision()), 
    //                                           &(ovl2.subdivision()) );
    //use_delete_notf = true;

    //int c_sweep_t;
    //c_sweep_t = clock();
    ovl->map_overlay(ovl1.subdivision(), ovl2.subdivision(), ovl_change_notf, arr_);
    //c_sweep_t = clock() - c_sweep_t;
    //std::cout<<"The time required by sweep line : "<< (double) c_sweep_t / (double) CLOCKS_PER_SEC<<std::endl;

    //sub_division1 = &ovl1;
    //sub_division2 = &ovl2;
  }
  
  Map_overlay (const Self &ovl1, const Self &ovl2, 
               Point_location_base *pl_ptr) : 
    arr_(pl_ptr), 
    sub_division1(&ovl1), sub_division2(&ovl2),
    ovl_change_notf(new Change_notification(&(ovl1.subdivision()), 
                                            &(ovl2.subdivision()) )), 
    ovl(new Map_ovl_sweep),
    use_delete_notf(true), use_delete_ovl(true)
  {
    //ovl = new Map_ovl_sweep;
    //use_delete_ovl = true;

    //ovl_change_notf = new Change_notification( &(ovl1.subdivision()), 
    //                                           &(ovl2.subdivision()) );
    //use_delete_notf = true;

    int c_sweep_t;
    c_sweep_t = clock();
    ovl->map_overlay(ovl1.subdivision(), ovl2.subdivision(), ovl_change_notf, arr_);
    c_sweep_t = clock() - c_sweep_t;
    std::cout<<"The time required by sweep line : "<< (double) c_sweep_t / (double) CLOCKS_PER_SEC<<std::endl;

    //sub_division1 = (Self *) &ovl1;
    //sub_division2 = (Self *) &ovl2;
  }

  Map_overlay (const Self &ovl1, 
               const Self &ovl2, 
               Change_notification* pmwx_change_notf) 
    : ovl_change_notf(pmwx_change_notf)
  {
    ovl = new Map_ovl_sweep;
    use_delete_ovl = true;
    use_delete_notf = false;

    ovl->map_overlay(ovl1.subdivision(), ovl2.subdivision(), ovl_change_notf, arr_);

    sub_division1 = (Self *) &ovl1;
    sub_division2 = (Self *) &ovl2;
  }

  Map_overlay (const Self &ovl1, const Self &ovl2, Map_ovl_base *ovl_ptr) 
    :  ovl(ovl_ptr) 
  {
    ovl_change_notf = new Change_notification( &(ovl1.subdivision()), 
                                                           &(ovl2.subdivision()) );
    use_delete_notf = true;
    use_delete_ovl = false;

    ovl->map_overlay(ovl1.subdivision(), ovl2.subdivision(), ovl_change_notf, arr_);
    
    sub_division1 = (Self *) &ovl1;
    sub_division2 = (Self *) &ovl2;
  }
  
  Map_overlay (const Self &ovl1, 
               const Self &ovl2, 
               Change_notification* pmwx_change_notf, 
               Map_ovl_base *ovl_ptr) 
    :  ovl_change_notf(pmwx_change_notf) , ovl(ovl_ptr) 
  {
    use_delete_notf = false;
    use_delete_ovl = false;

    ovl->map_overlay(ovl1.subdivision(), ovl2.subdivision(), ovl_change_notf, arr_);
    
    sub_division1 = (Self *) &ovl1;
    sub_division2 = (Self *) &ovl2;
  }
  
  virtual ~Map_overlay() 
  {
    if (use_delete_notf) delete ovl_change_notf;
    if (use_delete_ovl)  delete ovl;
  }

  // ------------------- Add a copy contr' 
  // -------------------- Add assignement operator.
  
  void	 delete_subtree() {
    arr_ = Map_overlay(arr_);
      
    for (Vertex_iterator v_iter = arr_.vertices_begin(); v_iter !=  arr_.vertices_end(); v_iter++)
      v_iter->reset();
    
    for (Halfedge_iterator h_iter = arr_.halfedges_begin(); h_iter !=  arr_.halfedges_end(); h_iter++)
      h_iter->reset();
    
    for (Face_iterator f_iter = arr_.faces_begin(); f_iter != arr_.faces_end(); f_iter++)
      f_iter->reset();
  }
  
  const Arrangement& subdivision() const 
  { return arr_; }
 
  const Self*  first_creator () const
  { return sub_division1; }
  
  const Self*  second_creator () const
  { return sub_division2; }

private:
  // only temporary.
  //typedef typename Arrangement::PMWXChangeNotification      PMWXChangeNotification_Empty;

  // copy a1 to a2.
  void  copy_arr(const Arrangement &a1, Arrangement &a2 /*, Map_overlay_ChangeNotification* ovl_change_notf*/)
  {
    typedef typename Arrangement::Traits    Traits;
    typedef typename Traits::X_curve        X_curve;
    
    cout<<"copy arrangment\n";

    for (Halfedge_const_iterator halfedge_iter = a1.halfedges_begin(); halfedge_iter != a1.halfedges_end(); halfedge_iter++, halfedge_iter++){
      //ovl_change_notf.set_curve_attributes(halfedge_iter->curve(), halfedge_iter, halfedge_iter->twin());
      //a2.insert(halfedge_iter->curve(), ovl_change_notf);
            
      a2.insert(halfedge_iter->curve());
    }

    // copy all attributes of halfedges and faces.
    Vertex_const_iterator v1_iter =  a1.vertices_begin();
    Vertex_iterator v2_iter = a2.vertices_begin();

    /*for ( ; v1_iter != a1.vertices_end() && v2_iter != a2.vertices_end(); v1_iter++, v2_iter++){
      v2_iter->set_bop(v1_iter->bop());
      }*/

    Halfedge_const_iterator h1_iter =  a1.halfedges_begin();
    Halfedge_iterator h2_iter = a2.halfedges_begin();

    for ( ; h1_iter != a1.halfedges_end() && h2_iter != a2.halfedges_end(); h1_iter++, h2_iter++){
      
      if (ovl_change_notf->get_first_halfedge_above(h1_iter) != h1_iter)
        ovl_change_notf->set_first_halfedge_above(h2_iter, ovl_change_notf->get_first_halfedge_above(h1_iter));

      if (ovl_change_notf->get_second_halfedge_above(h1_iter) != h1_iter)
        ovl_change_notf->set_second_halfedge_above(h2_iter, ovl_change_notf->get_second_halfedge_above(h1_iter));

      if (ovl_change_notf->get_first_face_above(h1_iter) != h1_iter->face())
        ovl_change_notf->set_first_face_above(h2_iter, ovl_change_notf->get_first_face_above(h1_iter));

      if (ovl_change_notf->get_second_face_above(h1_iter) != h1_iter->face())
        ovl_change_notf->set_second_face_above(h2_iter, ovl_change_notf->get_second_face_above(h1_iter));

      //h2_iter->set_bop(h1_iter->bop());  // defined only for bop.
    }

    Face_const_iterator f1_iter =  a1.faces_begin();
    Face_iterator f2_iter = a2.faces_begin();
    
    for ( ; f1_iter != a1.faces_end() && f2_iter != a2.faces_end(); f1_iter++, f2_iter++){

      if (ovl_change_notf->get_first_face_above(f1_iter) != f1_iter)
        ovl_change_notf->set_first_face_above(f2_iter, ovl_change_notf->get_first_face_above(f1_iter));

      if (ovl_change_notf->get_second_face_above(f1_iter) != f1_iter)
        ovl_change_notf->set_second_face_above(f2_iter, ovl_change_notf->get_second_face_above(f1_iter));
      
      //f2_iter->set_bop(f1_iter->bop());
    }
  }
  
  
  Arrangement                       arr_;
  const Self                        *sub_division1, *sub_division2;
  Change_notification               *ovl_change_notf;
  Map_ovl_base                      *ovl;
  bool use_delete_notf;
  bool use_delete_ovl;
};

CGAL_END_NAMESPACE

#endif







