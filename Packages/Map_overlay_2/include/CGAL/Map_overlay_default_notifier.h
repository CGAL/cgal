#ifndef CGAL_MAP_OVERLAY_DEFAULT_NOTIFIER_H
#define CGAL_MAP_OVERLAY_DEFAULT_NOTIFIER_H   

//#ifndef CGAL_MAP_OVERLAY_NAIVE_H
//#include <CGAL/Map_overlay_naive.h>
//#endif

CGAL_BEGIN_NAMESPACE         

template <class Arrangement_>
class Map_overlay_default_notifier : public Arrangement_::Change_notification
{
public:
  typedef Arrangement_                                                     Arrangement;
  typedef typename Arrangement::Vertex                                     Vertex;
  typedef typename Arrangement::Face                                       Face;
  typedef typename Arrangement::Halfedge                                   Halfedge;
  typedef typename Arrangement::Vertex_handle                              Vertex_handle;
  typedef typename Arrangement::Halfedge_handle                            Halfedge_handle;
  typedef typename Arrangement::Face_handle                                Face_handle;
  typedef typename Arrangement::Vertex_const_handle                        Vertex_const_handle;
  typedef typename Arrangement::Halfedge_const_handle                      Halfedge_const_handle;
  typedef typename Arrangement::Face_const_handle                          Face_const_handle;
  typedef typename Arrangement::Vertex_iterator                            Vertex_iterator;
  typedef typename Arrangement::Vertex_const_iterator                      Vertex_const_iterator;
  typedef typename Arrangement::Halfedge_iterator                          Halfedge_iterator;
  typedef typename Arrangement::Halfedge_const_iterator                    Halfedge_const_iterator;
  typedef typename Arrangement::Face_iterator                              Face_iterator;
  typedef typename Arrangement::Face_const_iterator                        Face_const_iterator;
  typedef typename Arrangement::Ccb_halfedge_circulator                    Ccb_halfedge_circulator;
  typedef typename Arrangement::Ccb_halfedge_const_circulator              Ccb_halfedge_const_circulator;
  typedef typename Arrangement::Holes_iterator                             Holes_iterator;
  typedef typename Arrangement::Holes_const_iterator                       Holes_const_iterator;
  typedef typename Arrangement::Locate_type                                Locate_type;
  typedef typename Arrangement::Traits_wrap                                Traits_wrap;
  typedef typename Arrangement::Change_notification                        Change_notification; 
  
  typedef typename Arrangement::Traits                                     Traits;
  typedef typename Traits::Point                                           Point;
  typedef typename Traits::X_curve                                         X_curve;

  typedef typename Arrangement::Planar_map                                PM;
  //typedef Arrangement_                                                   PM;
  typedef typename PM::Vertex_handle                                      Pm_vertex_handle;
  typedef typename PM::Halfedge_handle                                    Pm_halfedge_handle;
  typedef typename PM::Face_handle                                        Pm_face_handle;
  typedef typename PM::Vertex_const_handle                                Pm_vertex_const_handle;
  typedef typename PM::Halfedge_const_handle                              Pm_halfedge_const_handle;
  typedef typename PM::Face_const_handle                                  Pm_face_const_handle;

  typedef Map_overlay_default_notifier<Arrangement>        Self;
  typedef const Arrangement*                               Arr_const_pointer;
  
  Map_overlay_default_notifier() : arr1(0), arr2(0) {}

  Map_overlay_default_notifier (Arr_const_pointer sub_division1, 
                                Arr_const_pointer sub_division2) 
    : arr1(sub_division1), arr2(sub_division2) {}
  
  Map_overlay_default_notifier (const Self& notf) 
    : arr1(notf.arr1), arr2(notf.arr2) {}
  
  //Map_overlay_default_notifier (const Self* notf) 
  //  : arr1(notf->get_sub_division1()), arr2(notf->get_sub_division2()) {}
  
  virtual ~Map_overlay_default_notifier() {}

  void add_edge(const X_curve& cv, 
                Pm_halfedge_handle e, 
                bool left_to_right, 
                bool overlap = false)
  { 
    //std::cout<<"in add_edge" << std::endl;
    if ((CGAL::compare_lexicographically_xy(e->source()->point(), 
                                            e->target()->point()) == CGAL::SMALLER && 
         CGAL::compare_lexicographically_xy(orig_halfedge1->source()->point(), 
                                            orig_halfedge1->target()->point()) == CGAL::SMALLER) || 
        (CGAL::compare_lexicographically_xy(e->source()->point(), 
                                            e->target()->point()) == CGAL::LARGER && 
         CGAL::compare_lexicographically_xy(orig_halfedge1->source()->point(), 
                                            orig_halfedge1->target()->point()) == CGAL::LARGER)){
      if (first_halfedge){
        set_first_halfedge_above(Halfedge_handle(e), orig_halfedge1);
        set_first_halfedge_above(Halfedge_handle(e->twin()), orig_halfedge1->twin());

        //set_first_vertex_above(Vertex_handle(e->source()), orig_halfedge1->source());
        //set_first_vertex_above(Vertex_handle(e->target()), orig_halfedge1->target());
        
        set_first_halfedge_above(Vertex_handle(e->source()), orig_halfedge1);
        set_first_halfedge_above(Vertex_handle(e->target()), orig_halfedge1);
        
        //e->twin()->set_first_halfedge_above(orig_halfedge1->twin().operator->());
      }
      else{ 
        set_second_halfedge_above(Halfedge_handle(e), orig_halfedge1);
        set_second_halfedge_above(Halfedge_handle(e->twin()), orig_halfedge1->twin());

        //set_second_vertex_above(Vertex_handle(e->source()), orig_halfedge1->source());
        //set_second_vertex_above(Vertex_handle(e->target()), orig_halfedge1->target());
 
        set_second_halfedge_above(Vertex_handle(e->source()), orig_halfedge1);
        set_second_halfedge_above(Vertex_handle(e->target()), orig_halfedge1);

        //e->set_second_halfedge_above(orig_halfedge1.operator->());
        //e->twin()->set_second_halfedge_above(orig_halfedge1->twin().operator->());
      }
    }
    
    else{
      //if ((CGAL::compare_lexicographically_xy(e->source()->point(), e->target()->point()) == CGAL::SMALLER && 
      //              CGAL::compare_lexicographically_xy(orig_halfedge2->source()->point(), 
      //                                          orig_halfedge2->target()->point()) == CGAL::SMALLER)
      //      || (CGAL::compare_lexicographically_xy(e->source()->point(), 
      //                                             e->target()->point()) == CGAL::LARGER && 
      //          CGAL::compare_lexicographically_xy(orig_halfedge2->source()->point(), 
      //                                             orig_halfedge2->target()->point()) == CGAL::LARGER))
      if (first_halfedge){
        set_first_halfedge_above(Halfedge_handle(e), orig_halfedge2);
        set_first_halfedge_above(Halfedge_handle(e->twin()), orig_halfedge2->twin());
        
        //set_first_vertex_above(Vertex_handle(e->source()), orig_halfedge2->source());
        //set_first_vertex_above(Vertex_handle(e->target()), orig_halfedge2->target());
 
        set_first_halfedge_above(Vertex_handle(e->source()), orig_halfedge2);
        set_first_halfedge_above(Vertex_handle(e->target()), orig_halfedge2);

        //e->set_first_halfedge_above(orig_halfedge2.operator->());
        //e->twin()->set_first_halfedge_above(orig_halfedge2->twin().operator->());
      }
      else{
        set_second_halfedge_above(Halfedge_handle(e), orig_halfedge2);
        set_second_halfedge_above(Halfedge_handle(e->twin()), orig_halfedge2->twin());

        //set_second_vertex_above(Vertex_handle(e->source()), orig_halfedge2->source());
        //set_second_vertex_above(Vertex_handle(e->target()), orig_halfedge2->target());

        set_second_halfedge_above(Vertex_handle(e->source()), orig_halfedge2);
        set_second_halfedge_above(Vertex_handle(e->target()), orig_halfedge2);

        //e->set_second_halfedge_above(orig_halfedge2.operator->());
        //e->twin()->set_second_halfedge_above(orig_halfedge2->twin().operator->());
      }
    }

    // now making point location to update vertex pointers only if neccesary.
    if (get_first_vertex_above(e->source()) == e->source()){   // if the first vertex above the source of e is NULL.
      Point p = e->source()->point();
      Locate_type lt;
      Halfedge_const_handle h = arr1->locate(p, lt);
      if (lt == Arrangement::VERTEX){
        if (e->source()->point() ==  h->source()->point())
          set_first_vertex_above(Vertex_handle(e->source()), h->source());
        else
          set_first_vertex_above(Vertex_handle(e->source()), h->target());
      }
    }
    if (get_second_vertex_above(e->source()) == e->source()){   // if the second vertex above the source of e is NULL.
      Point p = e->source()->point();
      Locate_type lt;
      Halfedge_const_handle h = arr2->locate(p, lt);
      if (lt == Arrangement::VERTEX){
        if (e->source()->point() ==  h->source()->point())
          set_second_vertex_above(Vertex_handle(e->source()), h->source());
        else
          set_second_vertex_above(Vertex_handle(e->source()), h->target());
      }
    }
    
    if (get_first_vertex_above(e->target()) == e->target()){    // if the first vertex above the target of e is NULL.
      Point p = e->target()->point();
      Locate_type lt;
      Halfedge_const_handle h = arr1->locate(p, lt);
      if (lt == Arrangement::VERTEX){
        if (e->target()->point() ==  h->source()->point())
          set_first_vertex_above(Vertex_handle(e->target()), h->source());
        else
          set_first_vertex_above(Vertex_handle(e->target()), h->target());
      }
    }
    if (get_second_vertex_above(e->target()) == e->target()){   // if the second vertex above the target of e is NULL.
      Point p = e->target()->point();
      Locate_type lt;
      Halfedge_const_handle h = arr2->locate(p, lt);
      if (lt == Arrangement::VERTEX){
        if (e->target()->point() ==  h->source()->point())
          set_second_vertex_above(Vertex_handle(e->target()), h->source());
        else
          set_second_vertex_above(Vertex_handle(e->target()), h->target());
      }
    }

    // now making another test for the face to the new_edge side, 
    // in order to take care of overlapping halfedges.
    if ( !(e->face()->is_unbounded()) ){
      Ccb_halfedge_circulator  ccb_cir = Halfedge_handle(e)->face()->outer_ccb();
      
      do{
        if (ccb_cir->get_first_halfedge_above() != NULL){
          //std::cout<<"overlapping bogi"<<std::endl;
          set_first_face_above ( Halfedge_handle(e)->face(), 
                                 get_first_halfedge_above(ccb_cir)->face());
        }
        
        if (ccb_cir->get_second_halfedge_above() != NULL){ 
          set_second_face_above ( Halfedge_handle(e)->face(), 
                                  get_second_halfedge_above(ccb_cir)->face());
        }
        ++ccb_cir;
      } while (ccb_cir !=  Halfedge_handle(e)->face()->outer_ccb());
    }
    
    if ( !(e->twin()->face()->is_unbounded()) ){
      Ccb_halfedge_circulator  ccb_cir = Halfedge_handle(e->twin())->face()->outer_ccb();
      
      do{
        if (ccb_cir->get_first_halfedge_above() != NULL){
          //std::cout<<"overlapping twin bogi"<<std::endl;
          set_first_face_above ( Halfedge_handle(e->twin())->face(), 
                                 get_first_halfedge_above(ccb_cir)->face());
        }
        
        if (ccb_cir->get_second_halfedge_above() != NULL){ 
          set_second_face_above ( Halfedge_handle(e->twin())->face(), 
                                  get_second_halfedge_above(ccb_cir)->face());
        }
        ++ccb_cir;
      } while (ccb_cir !=  Halfedge_handle(e->twin())->face()->outer_ccb());
    }
  }
  
  void split_edge(Pm_halfedge_handle orig_edge, 
                  Pm_halfedge_handle new_edge, 
                  const X_curve& c1, 
                  const X_curve& c2)
  {
    // update half edges above new_edge and its twin.
    new_edge->set_first_halfedge_above(orig_edge->get_first_halfedge_above());
    new_edge->set_second_halfedge_above(orig_edge->get_second_halfedge_above());
    new_edge->twin()->set_first_halfedge_above(orig_edge->twin()->get_first_halfedge_above());
    new_edge->twin()->set_second_halfedge_above(orig_edge->twin()->get_second_halfedge_above());
    
    // upadate halfedge above the edge points of new_edge.
    if (get_first_halfedge_above(orig_edge) != orig_edge){
      set_first_halfedge_above(Vertex_handle(new_edge->source()), get_first_halfedge_above(orig_edge));
      set_first_halfedge_above(Vertex_handle(new_edge->target()), get_first_halfedge_above(orig_edge));
    }
    if (get_second_halfedge_above(orig_edge) != orig_edge){
      set_second_halfedge_above(Vertex_handle(new_edge->source()), get_second_halfedge_above(orig_edge));
      set_second_halfedge_above(Vertex_handle(new_edge->target()), get_second_halfedge_above(orig_edge));
    }
    
    // now making a nother test for the face to the new_edge side, in order to take care of overlapping halfedges.
    if (new_edge->face()->is_unbounded())
      return;
    
    Ccb_halfedge_circulator  ccb_cir = Halfedge_handle(new_edge)->face()->outer_ccb();

    do{
      if (ccb_cir->get_first_halfedge_above() != NULL){
        //std::cout<<"overlapping bogi"<<std::endl;
        set_first_face_above ( Halfedge_handle(new_edge)->face(), 
                               get_first_halfedge_above(ccb_cir)->face());
      }
      
      if (ccb_cir->get_second_halfedge_above() != NULL){ 
        set_second_face_above ( Halfedge_handle(new_edge)->face(), 
                                get_second_halfedge_above(ccb_cir)->face());
      }
      ++ccb_cir;
    } while (ccb_cir !=  Halfedge_handle(new_edge)->face()->outer_ccb());
  }
  
  void update_splited_face (Face_handle new_face)
  {
    new_face->reset();

    // taking care of the unbounded face.
    if (new_face->is_unbounded()){
      set_first_face_above(new_face, arr1->unbounded_face());
      set_second_face_above(new_face, arr2->unbounded_face());
      return;
    }

    Ccb_halfedge_circulator  ccb_cir = new_face->outer_ccb();
    
    //std::cout<<"updating face"<<std::endl;

    do{
      if (ccb_cir->get_first_halfedge_above() != NULL){ 
        set_first_face_above (new_face, get_first_halfedge_above(ccb_cir)->face());
        //std::cout<<ccb_cir->curve()<<" had blue face above "<<std::endl;
        //new_face->set_first_face_above( ((Halfedge* ) ccb_cir->get_first_halfedge_above())->face().operator->());
      }
      
      if (ccb_cir->get_second_halfedge_above() != NULL){ 
        set_second_face_above (new_face, get_second_halfedge_above(ccb_cir)->face());
        //std::cout<<ccb_cir->curve()<<" had red face above "<<std::endl;
        //new_face->set_second_face_above(((Halfedge* ) ccb_cir->get_second_halfedge_above())->face().operator->());
      }
      ++ccb_cir;
    } while (ccb_cir != new_face->outer_ccb());

    // now updating the faces above the halfedges and vertices of the new face.
    /*do{
      if (new_face->get_first_face_above() != NULL){ 
        set_first_face_above (ccb_cir, get_first_face_above(new_face));
        set_first_face_above (ccb_cir->source(), get_first_face_above(new_face));
      }
      
      if (new_face->get_second_face_above() != NULL){ 
        set_second_face_above (ccb_cir, get_second_face_above(new_face));
        set_second_face_above (ccb_cir->source(), get_second_face_above(new_face));
      }
      ++ccb_cir;
      } while (ccb_cir != new_face->outer_ccb());*/
    
    // making point location to the original subdivision 
    // (if new_face is a hole of at least one of them).
    if (new_face->get_first_face_above() == NULL){
      ccb_cir != new_face->outer_ccb();
      do{
        Point p = ccb_cir->source()->point();
        //std::cout<<"internal point is : "<< p <<std::endl;
        
        Locate_type lt;
        Halfedge_const_handle e = arr1->locate(p, lt);
        if (lt == Arrangement::FACE || lt == Arrangement::UNBOUNDED_FACE)
          {
            //cout<<"internal point in FACE is : "<< p <<endl;
            //if (e->face()->is_unbounded())
            //  cout<<"Unbounded"<<endl;
            //else
            //  cout<<"Bounded"<<endl;

            set_first_face_above(new_face, e->face());
            set_first_face_above(ccb_cir, e->face());
            set_first_face_above(ccb_cir->source(), e->face());  // also set the face above the vertex.
            break;
          }
        ++ccb_cir;
      } while (ccb_cir != new_face->outer_ccb());
    }
    
    if (new_face->get_second_face_above() == NULL){
      ccb_cir != new_face->outer_ccb();
      do{
        Point p = ccb_cir->source()->point();
        //cout<<"internal point is : "<< p <<"\n";
        
        Locate_type lt;
        Halfedge_const_handle e = arr2->locate(p, lt);
        if (lt == Arrangement::FACE || lt == Arrangement::UNBOUNDED_FACE)
          {
            //cout<<"internal point in FACE is : "<< p <<endl;
            //if (e->face()->is_unbounded())
            //  cout<<"Unbounded"<<endl;
            //else
            //  cout<<"Bounded"<<endl;
            
            set_second_face_above(new_face, e->face());
            set_second_face_above(ccb_cir, e->face());
            set_second_face_above(ccb_cir->source(), e->face());  // also set the face above the vertex.
            break;
          }
        ++ccb_cir;
      } while (ccb_cir != new_face->outer_ccb());
    }

    // now updating the faces above the halfedges and vertices of the new face.
    do{
      if (new_face->get_first_face_above() != NULL){ 
        set_first_face_above (ccb_cir, get_first_face_above(new_face));
        set_first_face_above (ccb_cir->source(), get_first_face_above(new_face));
      }
      
      if (new_face->get_second_face_above() != NULL){ 
        set_second_face_above (ccb_cir, get_second_face_above(new_face));
        set_second_face_above (ccb_cir->source(), get_second_face_above(new_face));
      }
      ++ccb_cir;
    } while (ccb_cir != new_face->outer_ccb());
  }
  
  void split_face(Pm_face_handle orig_face, Pm_face_handle new_face)
  {

    //std::cout<<"is split_face"<<std::endl;

    //if (orig_face == new_face)
    //  return; 

    update_splited_face (Face_handle(orig_face));
    
    update_splited_face (Face_handle(new_face));
  }
  
  
  void add_hole(Pm_face_handle in_face, Pm_halfedge_handle new_hole) {}
  
  
  const X_curve &edge_support_curve(Halfedge_handle edge)
  {
    return edge->curve();
  }
  
  bool have_support_curve()
  {
    return false;
  }
  
  /***************************************** new functions **********************************************************/
  void set_curve_attributes(const X_curve& cv, 
                            Halfedge_const_handle orig_halfedge1_, 
                            //Halfedge_const_handle orig_halfedge2_, 
                            bool first_halfedge_)
  {
    orig_halfedge1 = orig_halfedge1_;
    orig_halfedge2 = orig_halfedge1_->twin();
    first_halfedge = first_halfedge_;
  }
  
  Arr_const_pointer first_subdivision () const { return arr1;}
  
  Arr_const_pointer second_subdivision () const { return arr2;}
  
  //-----------------------------------------  handle wrappering.
  // setting the vertex above.
  void  set_first_vertex_above(Vertex_handle v, Vertex_handle v_above) {
    v->set_first_vertex_above(v_above.operator->());
  }

  void  set_second_vertex_above(Vertex_handle v, Vertex_handle v_above) {
    v->set_second_vertex_above(v_above.operator->());
  }

  void  set_first_vertex_above(Vertex_handle v, Vertex_const_handle v_above) {
    v->set_first_vertex_above(v_above.operator->());
  }

  void  set_second_vertex_above(Vertex_handle v, Vertex_const_handle v_above) {
    v->set_second_vertex_above(v_above.operator->());
  }

  // setting the halfedge above.
  void  set_first_halfedge_above(Vertex_handle v, Halfedge_handle h_above) {
    v->set_first_halfedge_above(h_above.operator->());
  }

  void  set_second_halfedge_above(Vertex_handle v, Halfedge_handle h_above) {
    v->set_second_halfedge_above(h_above.operator->());
  }

  void  set_first_halfedge_above(Vertex_handle v, Halfedge_const_handle h_above) {
    v->set_first_halfedge_above(h_above.operator->());
  }
  
  void  set_second_halfedge_above(Vertex_handle v, Halfedge_const_handle h_above) {
    v->set_second_halfedge_above(h_above.operator->());
  }

  void  set_first_halfedge_above(Halfedge_handle h, Halfedge_handle h_above) {
    h->set_first_halfedge_above(h_above.operator->());
  }

  void  set_second_halfedge_above(Halfedge_handle h, Halfedge_handle h_above) {
    h->set_second_halfedge_above(h_above.operator->());
  }

  void  set_first_halfedge_above(Halfedge_handle h, Halfedge_const_handle h_above) {
    h->set_first_halfedge_above(h_above.operator->());
  }
  
  void  set_second_halfedge_above(Halfedge_handle h, Halfedge_const_handle h_above) {
    h->set_second_halfedge_above(h_above.operator->());
  }
  
  // setting the faces above.
  void  set_first_face_above(Vertex_handle v, Face_handle f_above) {
    v->set_first_face_above(f_above.operator->());
  }

  void  set_second_face_above(Vertex_handle v, Face_handle f_above) {
    v->set_second_face_above(f_above.operator->());
  }
  
  void  set_first_face_above(Vertex_handle v, Face_const_handle f_above) {
    v->set_first_face_above(f_above.operator->());
  }
  
  void  set_second_face_above(Vertex_handle v, Face_const_handle f_above) {
    v->set_second_face_above(f_above.operator->());
  }

  void  set_first_face_above(Halfedge_handle h, Face_handle f_above) {
    h->set_first_face_above(f_above.operator->());
  }

  void  set_second_face_above(Halfedge_handle h, Face_handle f_above) {
    h->set_second_face_above(f_above.operator->());
  }

  void  set_first_face_above(Halfedge_handle h, Face_const_handle f_above) {
    h->set_first_face_above(f_above.operator->());
  }
  
  void  set_second_face_above(Halfedge_handle h, Face_const_handle f_above) {
    h->set_second_face_above(f_above.operator->());
  }
  
  void  set_first_face_above(Face_handle f, Face_handle f_above) {
    f->set_first_face_above(f_above.operator->());
  }

  void  set_second_face_above(Face_handle f, Face_handle f_above) {
    f->set_second_face_above(f_above.operator->());
  }

  void  set_first_face_above(Face_handle f, Face_const_handle f_above) {
    f->set_first_face_above(f_above.operator->());
  }
  
  void  set_second_face_above(Face_handle f, Face_const_handle f_above) {
    f->set_second_face_above(f_above.operator->());
  }

  // getting the vertex above.
  Vertex_handle get_first_vertex_above(Vertex_handle v) {
    Vertex* vp = (Vertex*) v->get_first_vertex_above() ;
    
    if (vp){
      Pm_vertex_handle tmp_v = Pm_vertex_handle(vp);
      return Vertex_handle(tmp_v);
    }
    else
      return v;
  }
    
  Vertex_handle get_second_vertex_above(Vertex_handle v) {
    Vertex* vp = (Vertex*) v->get_second_vertex_above() ;
    
    if (vp){
      Pm_vertex_handle tmp_v = Pm_vertex_handle(vp);
      return Vertex_handle(tmp_v);
    }
    else
      return v;
  }
  
  Vertex_const_handle get_first_vertex_above(Vertex_const_handle v) {
    Vertex* vp = (Vertex*) v->get_first_vertex_above() ;
    
    if (vp){
      Pm_vertex_const_handle tmp_v = Pm_vertex_const_handle(vp);
      return Vertex_const_handle(tmp_v);
    }
    else
      return v;
  }
  
  Vertex_const_handle get_second_vertex_above(Vertex_const_handle v) {
    Vertex* vp = (Vertex*) v->get_second_vertex_above() ;
    
    if (vp){
      Pm_vertex_const_handle tmp_v = Pm_vertex_const_handle(vp);
      return Vertex_const_handle(tmp_v);
    }
    else
      return v;
  }
  
  // getting the halfedge above.
  Halfedge_handle get_first_halfedge_above(Vertex_handle v) {
    Halfedge* hp = (Halfedge*) v->get_first_halfedge_above() ;
    
    if (hp){
      Pm_halfedge_handle tmp_h = Pm_halfedge_handle(hp);
      return Halfedge_handle(tmp_h);
    }
    else
      return v->incident_halfedges();
  }
   
  Halfedge_handle get_second_halfedge_above(Vertex_handle v) {
    Halfedge* hp = (Halfedge*) v->get_second_halfedge_above() ;

    if (hp){
      Pm_halfedge_handle tmp_h = Pm_halfedge_handle(hp);
      return Halfedge_handle(tmp_h);
    }
    else
      return v->incident_halfedges();
  }
  
  Halfedge_const_handle get_first_halfedge_above(Vertex_const_handle v) {
    Halfedge* hp = (Halfedge*) v->get_first_halfedge_above() ;
    
    if (hp){
      Pm_halfedge_const_handle tmp_h = Pm_halfedge_const_handle(hp);
      return Halfedge_const_handle(tmp_h);
    }
    else
      return v->incident_halfedges();
  }
   
  Halfedge_const_handle get_second_halfedge_above(Vertex_const_handle v) {
    Halfedge* hp = (Halfedge*) v->get_second_halfedge_above() ;
    
    if (hp){
      Pm_halfedge_const_handle tmp_h = Pm_halfedge_const_handle(hp);
      return Halfedge_const_handle(tmp_h);
    }
    else
      return v->incident_halfedges();
  }

  Halfedge_handle get_first_halfedge_above(Halfedge_handle h) {
    Halfedge* hp = (Halfedge*) h->get_first_halfedge_above() ;
    
    if (hp){
      Pm_halfedge_handle tmp_h = Pm_halfedge_handle(hp);
      return Halfedge_handle(tmp_h);
    }
    else
      return h;
  }
   
  Halfedge_handle get_second_halfedge_above(Halfedge_handle h) {
    Halfedge* hp = (Halfedge*) h->get_second_halfedge_above() ;

    if (hp){
      Pm_halfedge_handle tmp_h = Pm_halfedge_handle(hp);
      return Halfedge_handle(tmp_h);
    }
    else
      return h;
  }

  Halfedge_const_handle get_first_halfedge_above(Halfedge_const_handle h) {
    Halfedge* hp = (Halfedge*) h->get_first_halfedge_above() ;
    
    if (hp){
      Pm_halfedge_const_handle tmp_h = Pm_halfedge_const_handle(hp);
      return Halfedge_const_handle(tmp_h);
    }
    else
      return h;
  }
   
  Halfedge_const_handle get_second_halfedge_above(Halfedge_const_handle h) {
    Halfedge* hp = (Halfedge*) h->get_second_halfedge_above() ;
    
    if (hp){
      Pm_halfedge_const_handle tmp_h = Pm_halfedge_const_handle(hp);
      return Halfedge_const_handle(tmp_h);
    }
    else
      return h;
  }

  // getting the face above.

  Face_handle get_first_face_above(Vertex_handle v) {
    Face* fp = (Face*) v->get_first_face_above() ;
    
    if (fp){
      Pm_face_handle tmp_f = Pm_face_handle(fp);
      return Face_handle(tmp_f);
    }
    else
      return v->incident_halfedges()->face();
  }
  
  Face_handle get_second_face_above(Vertex_handle v) {
    Face* fp = (Face*) v->get_second_face_above() ;
    
    if (fp){
      Pm_face_handle tmp_f = Pm_face_handle(fp);
      return Face_handle(tmp_f);
    }
    else
      return v->incident_halfedges()->face();
  }

  Face_const_handle get_first_face_above(Vertex_const_handle v) {
    Face* fp=(Face*) v->get_first_face_above() ;
    
    if (fp){
      Pm_face_const_handle tmp_f = Pm_face_const_handle(fp);
      return Face_const_handle(tmp_f);
    }
    else
      return v->incident_halfedges()->face();
  }
  
  Face_const_handle get_second_face_above(Vertex_const_handle v) {
    Face* fp=(Face*) v->get_second_face_above() ;
    
    if (fp){
      Pm_face_const_handle tmp_f = Pm_face_const_handle(fp);
      return Face_const_handle(tmp_f);
    }
    else
      return v->incident_halfedges()->face();
  }

  Face_handle get_first_face_above(Halfedge_handle h) {
    Face* fp = (Face*) h->get_first_face_above() ;
    
    if (fp){
      Pm_face_handle tmp_f = Pm_face_handle(fp);
      return Face_handle(tmp_f);
    }
    else
      return h->face();
  }
  
  Face_handle get_second_face_above(Halfedge_handle h) {
    Face* fp = (Face*) h->get_second_face_above() ;
    
    if (fp){
      Pm_face_handle tmp_f = Pm_face_handle(fp);
      return Face_handle(tmp_f);
    }
    else
      return h->face();
  }

  Face_const_handle get_first_face_above(Halfedge_const_handle h) {
    Face* fp=(Face*) h->get_first_face_above() ;
    
    if (fp){
      Pm_face_const_handle tmp_f = Pm_face_const_handle(fp);
      return Face_const_handle(tmp_f);
    }
    else
      return h->face();
  }
  
  Face_const_handle get_second_face_above(Halfedge_const_handle h) {
    Face* fp=(Face*) h->get_second_face_above() ;
    
    if (fp){
      Pm_face_const_handle tmp_f = Pm_face_const_handle(fp);
      return Face_const_handle(tmp_f);
    }
    else
      return h->face();
  }

  Face_handle get_first_face_above(Face_handle f) {
    Face* fp = (Face*) f->get_first_face_above() ;
    
    if (fp){
      Pm_face_handle tmp_f = Pm_face_handle(fp);
      return Face_handle(tmp_f);
    }
    else
      return f;
  }
  
  Face_handle get_second_face_above(Face_handle f) {
    Face* fp = (Face*) f->get_second_face_above() ;
    
    if (fp){
      Pm_face_handle tmp_f = Pm_face_handle(fp);
      return Face_handle(tmp_f);
    }
    else
      return f;
  }
  
  Face_const_handle get_first_face_above(Face_const_handle f) {
    Face* fp = (Face*) f->get_first_face_above() ;
    
    if (fp){
      Pm_face_const_handle tmp_f = Pm_face_const_handle(fp);
      return Face_const_handle(tmp_f);
    }
    else
      return f;
  }
  
  Face_const_handle get_second_face_above(Face_const_handle f) {
    Face* fp=(Face*) f->get_second_face_above() ;
 
    if (fp){
      Pm_face_const_handle tmp_f = Pm_face_const_handle(fp);
      return Face_const_handle(tmp_f);
    }
    else
      return f;
  }
   
  //---------------------------- end handles wrappering.
  
private:
  bool                   first_halfedge;
  Halfedge_const_handle  orig_halfedge1, orig_halfedge2;
  Arr_const_pointer      arr1, arr2;
};

CGAL_END_NAMESPACE

#endif



