#ifndef CGAL_MAP_OVERLAY_H
#define CGAL_MAP_OVERLAY_H

#include <vector>
#include <list>
//#include <time.h>

#ifndef CGAL_PM_POINT_LOCATION_BASE_H
#include <CGAL/Pm_point_location_base.h>
#endif

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

template  <class Subdivision_, 
           class Change_notification_ = Map_overlay_default_notifier<Subdivision_> > 
class  Map_overlay {
public:  
  
  typedef Subdivision_                                          Subdivision;
  typedef Change_notification_                                  Change_notification;
  typedef typename Subdivision::Planar_map                      Planar_map;
  //typedef typename Subdivision::Pmwx                          Pmwx;
  typedef typename Subdivision::Vertex                          Vertex;
  typedef typename Subdivision::Face                            Face;
  typedef typename Subdivision::Halfedge                        Halfedge;
  typedef typename Subdivision::Vertex_handle                   Vertex_handle;
  typedef typename Subdivision::Halfedge_handle                 Halfedge_handle;
  typedef typename Subdivision::Face_handle                     Face_handle;
  typedef typename Subdivision::Vertex_const_handle             Vertex_const_handle;
  typedef typename Subdivision::Halfedge_const_handle           Halfedge_const_handle;
  typedef typename Subdivision::Face_const_handle               Face_const_handle;
  typedef typename Subdivision::Vertex_iterator                 Vertex_iterator;
  typedef typename Subdivision::Halfedge_iterator               Halfedge_iterator;
  typedef typename Subdivision::Face_iterator                   Face_iterator;
  typedef typename Subdivision::Vertex_const_iterator           Vertex_const_iterator;
  typedef typename Subdivision::Halfedge_const_iterator         Halfedge_const_iterator;
  typedef typename Subdivision::Face_const_iterator             Face_const_iterator;
  typedef typename Subdivision::Ccb_halfedge_circulator         Ccb_halfedge_circulator;
  typedef typename Subdivision::Ccb_halfedge_const_circulator   Ccb_halfedge_const_circulator;
  typedef typename Subdivision::Holes_iterator                  Holes_iterator;
  typedef typename Subdivision::Holes_const_iterator            Holes_const_iterator;
  typedef typename Subdivision::Locate_type                     Locate_type;
  //typedef typename Subdivision::Point_location_base             Point_location_base;

  typedef typename Subdivision::Traits                         Traits;
  typedef typename Traits::X_curve_2                           X_curve_2;  
  typedef typename Traits::Point_2                             Point_2;

  typedef Map_overlay_base<Subdivision, Change_notification>   Map_overlay_base;
  typedef Map_overlay_sweep<Subdivision, Change_notification>  Map_overlay_sweep;
  
  typedef Pm_point_location_base<Planar_map> Point_location_base;
  
private:
  typedef Map_overlay<Subdivision, Change_notification>       Self;
  
public:

  Map_overlay() : 
    arr_(new Subdivision()),
    sub_division1(0), sub_division2(0), 
    ovl_change_notf(new Change_notification),  
    ovl_alg(new Map_overlay_sweep), 
    use_delete_notf(true), 
    use_delete_ovl(true) {}

   Map_overlay(Point_location_base *pl_ptr) : 
     arr_(new Subdivision(pl_ptr)), 
     sub_division1(0), sub_division2(0), 
     ovl_change_notf(new Change_notification),  
     ovl_alg(new Map_overlay_sweep), 
     use_delete_notf(true), 
     use_delete_ovl(true) {}
  
  Map_overlay (const Subdivision &arr) : 
    arr_(new Subdivision(arr)), 
    sub_division1(0), sub_division2(0),  
    ovl_change_notf(new Change_notification),  
    ovl_alg(new Map_overlay_sweep), 
    use_delete_notf(true), 
    use_delete_ovl(true) {}
  
  Map_overlay (const Subdivision &arr, 
               Change_notification* pmwx_change_notf) : 
    arr_(new Subdivision(arr)), 
    sub_division1(0), sub_division2(0), 
    ovl_change_notf(pmwx_change_notf), 
    ovl_alg(new Map_overlay_sweep), 
    use_delete_notf(false), use_delete_ovl(true) {}
  
  Map_overlay (const Subdivision &arr, 
               Map_overlay_base *ovl_ptr) : 
    arr_(new Subdivision(arr)),  
    sub_division1(0), sub_division2(0), 
    ovl_change_notf(new Change_notification), 
    ovl_alg(ovl_ptr), 
    use_delete_notf(true), use_delete_ovl(false) {}

  Map_overlay (const Subdivision &arr,  
               Change_notification* pmwx_change_notf, 
               Map_overlay_base *ovl_ptr) : 
    arr_(new Subdivision(arr)), 
    sub_division1(0), sub_division2(0), 
    ovl_change_notf(pmwx_change_notf), ovl_alg(ovl_ptr),
    use_delete_notf(false), use_delete_ovl(false) {}

  Map_overlay (const Self &ovl1, const Self &ovl2) : 
    arr_(new Subdivision),
    sub_division1(&ovl1), sub_division2(&ovl2),
    ovl_change_notf(new Change_notification(&(ovl1.subdivision()), 
                                            &(ovl2.subdivision()) )), 
    ovl_alg(new Map_overlay_sweep), 
    use_delete_notf(true), use_delete_ovl(true)
  {
    ovl_alg->map_overlay(ovl1.subdivision(), ovl2.subdivision(), ovl_change_notf, *arr_);
  }
  
  Map_overlay (const Self &ovl1, const Self &ovl2, 
               Point_location_base *pl_ptr) : 
    arr_(new Subdivision(pl_ptr)), 
    sub_division1(&ovl1), sub_division2(&ovl2),
    ovl_change_notf(new Change_notification(&(ovl1.subdivision()), 
                                            &(ovl2.subdivision()) )), 
    ovl_alg(new Map_overlay_sweep),
    use_delete_notf(true), use_delete_ovl(true)
  {
    //int c_sweep_t;
    //c_sweep_t = clock();
    
    ovl_alg->map_overlay(ovl1.subdivision(), ovl2.subdivision(), ovl_change_notf, *arr_);
    
    //c_sweep_t = clock() - c_sweep_t;
    //std::cout<<"The time required by sweep line : "<< (double) c_sweep_t / (double) CLOCKS_PER_SEC<<std::endl;
  }

  Map_overlay (const Self &ovl1, 
               const Self &ovl2, 
               Change_notification* pmwx_change_notf) : 
    arr_(new Subdivision()),
    sub_division1(&ovl1), sub_division2(&ovl2),
    ovl_change_notf(pmwx_change_notf), ovl_alg(new Map_overlay_sweep),
    use_delete_notf(false), use_delete_ovl(true)
  {
    ovl_alg->map_overlay(ovl1.subdivision(), ovl2.subdivision(), ovl_change_notf, *arr_);
  }

   Map_overlay (const Self &ovl1, 
                const Self &ovl2, 
                Point_location_base* pl_ptr,
                Change_notification* pmwx_change_notf) :  
     arr_(new Subdivision(pl_ptr)), 
     sub_division1(&ovl1), sub_division2(&ovl2),
     ovl_change_notf(pmwx_change_notf), ovl_alg(new Map_overlay_sweep),
     use_delete_notf(false), use_delete_ovl(true)
  {
    ovl_alg->map_overlay(ovl1.subdivision(), ovl2.subdivision(), ovl_change_notf, *arr_);
  }

  
  Map_overlay (const Self &ovl1, 
               const Self &ovl2, 
               Map_overlay_base *ovl_ptr) :
    arr_(new Subdivision()),
    sub_division1(&ovl1), sub_division2(&ovl2),
    ovl_change_notf(new Change_notification(&(ovl1.subdivision()), 
                                            &(ovl2.subdivision()) )),
    ovl_alg(ovl_ptr), 
    use_delete_notf(true), use_delete_ovl(false)
  {
    ovl_alg->map_overlay(ovl1.subdivision(), ovl2.subdivision(), ovl_change_notf, *arr_);
  }
  
  Map_overlay (const Self &ovl1, 
               const Self &ovl2, 
               Point_location_base* pl_ptr, 
               Map_overlay_base *ovl_ptr) :
    arr_(new Subdivision(pl_ptr)),
    sub_division1(&ovl1), sub_division2(&ovl2),
    ovl_change_notf(new Change_notification(&(ovl1.subdivision()), 
                                            &(ovl2.subdivision()) )),
    ovl_alg(ovl_ptr), 
    use_delete_notf(true), use_delete_ovl(false)
  {
    ovl_alg->map_overlay(ovl1.subdivision(), ovl2.subdivision(), ovl_change_notf, *arr_);
  }
  
  
  Map_overlay (const Self &ovl1, 
               const Self &ovl2, 
               Change_notification* pmwx_change_notf, 
               Map_overlay_base *ovl_ptr) :
    arr_(new Subdivision()),
    sub_division1(&ovl1), sub_division2(&ovl2),
    ovl_change_notf(pmwx_change_notf), ovl_alg(ovl_ptr),
    use_delete_notf(false), use_delete_ovl(false)
  {
    ovl_alg->map_overlay(ovl1.subdivision(), ovl2.subdivision(), ovl_change_notf, *arr_);
  }

  Map_overlay (const Self &ovl1, 
               const Self &ovl2, 
               Point_location_base* pl_ptr,
               Change_notification* pmwx_change_notf, 
               Map_overlay_base *ovl_ptr) :
    arr_(new Subdivision(pl_ptr)),
    sub_division1(&ovl1), sub_division2(&ovl2),
    ovl_change_notf(pmwx_change_notf), ovl_alg(ovl_ptr),
    use_delete_notf(false), use_delete_ovl(false)
  {
    ovl_alg->map_overlay(ovl1.subdivision(), ovl2.subdivision(), ovl_change_notf, *arr_);
  }
  
  // ------------------- Copy contructor -----------------------------
  Map_overlay (const Self &ovl) : 
    arr_(new Subdivision(*(ovl.arr_))), 
    sub_division1(ovl.sub_division1), 
    sub_division2(ovl.sub_division2), 
    ovl_change_notf(new Change_notification(ovl.ovl_change_notf)),
    ovl_alg(ovl.ovl_alg),  // needs a copy constructor for Point location.
    use_delete_notf(true), use_delete_ovl(false) {}
  
  virtual ~Map_overlay() 
  {
    delete arr_;
    if (use_delete_notf) delete ovl_change_notf;
    if (use_delete_ovl)  delete ovl_alg;
  }
  
  // -------------------- Assignement operator --------------------------
  const Self& operator=(const Self& ovl)
  {
    *arr_ = *(ovl.arr_); 
    sub_division1 = ovl.sub_division1; 
    sub_division2 = ovl.sub_division2; 

    // The notifier and ovl algorithm remain their initial values obtained 
    // in the defualt constructor (or other constructor taken for *this).
    //ovl_change_notf(new Change_notification(ovl.ovl_change_notf));
    //ovl_alg = ovl.ovl_alg;
    //use_delete_notf=true; use_delete_ovl=false;

    return *this;
  }

  void	 delete_subtree() {
    (*arr_) = Map_overlay(*arr_);
      
    for (Vertex_iterator v_iter = arr_->vertices_begin(); 
         v_iter !=  arr_->vertices_end(); v_iter++)
      v_iter->reset();
    
    for (Halfedge_iterator h_iter = arr_->halfedges_begin(); 
         h_iter !=  arr_->halfedges_end(); h_iter++)
      h_iter->reset();
    
    for (Face_iterator f_iter = arr_->faces_begin(); 
         f_iter != arr_->faces_end(); f_iter++)
      f_iter->reset();
  }
  
  const Subdivision& subdivision() const 
  { return *arr_; }
 
  const Self*  first_creator () const
  { return sub_division1; }
  
  const Self*  second_creator () const
  { return sub_division2; }

  const Change_notification* change_notification() const 
  { return ovl_change_notf; }
  
  
private:
  // Saving the subdivision as a pointer is crucial.
  // Since the overlay has levels of creators, we make sure this way all levels are valid.
  // For example, it is possible to create a map overlay of two subdivisions (rather than
  // two map overlays) using the converter and no harm is done: the converter saves its 
  // corresponding map overlay as the creator. If the usage is not by pointers, a temporary 
  // value of map overlay may be constructed by the converter and may cause memory fligh.
  Subdivision                       *arr_;
  const Self                        *sub_division1, *sub_division2;
  Change_notification               *ovl_change_notf;
  Map_overlay_base                  *ovl_alg;
  bool use_delete_notf;
  bool use_delete_ovl;
};

CGAL_END_NAMESPACE

#endif







