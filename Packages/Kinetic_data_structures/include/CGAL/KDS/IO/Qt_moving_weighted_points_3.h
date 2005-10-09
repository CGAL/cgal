#ifndef CGAL_KDS_QT_MOVING_WEIGHTED_POINT_TABLE_3_H
#define CGAL_KDS_QT_MOVING_WEIGHTED_POINT_TABLE_3_H
#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/Regular_triangulation_instantaneous_traits_3.h>
#include <CGAL/KDS/Ref_counted.h>
#include <CGAL/KDS/Simulator_objects_listener.h>
#include <CGAL/KDS/IO/Coin_pointer.h>

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodekits/SoShapeKit.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoPointSet.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoText2.h>
#include <Inventor/nodekits/SoAppearanceKit.h>
//#include <CGAL/KDS/IO/Qt_gui_3.h>

#include <Inventor/events/SoEvent.h>
#include <Inventor/events/SoButtonEvent.h>
#include <Inventor/events/SoKeyboardEvent.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/sensors/SoOneShotSensor.h>

CGAL_KDS_BEGIN_NAMESPACE;

//! A graphical moving point set in 3D
/*!
  This class shows how to interact with the Coin and Qt based 3D gui. 

	How points are drawn can be controlled from the keyboard
	- L toggles drawing labels.
	- S has the points be drawn as spheres with radius equal to the square root of the weight.
	- P has the points be drawn as GL points with some (constant) size
	- T has tbe points be drawn as transparent spheres
*/
template <class Traits, class GUI> 
class Qt_moving_weighted_points_3: public Ref_counted< Qt_moving_weighted_points_3<Traits, GUI> > {
protected:
  typedef Qt_moving_weighted_points_3<Traits, GUI> This;
  typedef typename Traits::Instantaneous_kernel IK;
  typedef typename Traits::Moving_point_table MPT;
  typedef typename Traits::Simulator Simulator;
  typedef typename GUI::Listener Gui_listener;
  typedef typename Simulator::Listener Simulator_listener;
  typedef CGAL::KDS::Simulator_objects_listener<Simulator_listener, This> Siml;
  friend class CGAL::KDS::Simulator_objects_listener<Simulator_listener, This>;
  class Guil;

  typedef typename MPT::Data::Coordinate::NT NT;
  
  class Table_listener: public MPT::Listener{
    typedef typename MPT::Listener P;
  public:
    Table_listener(typename MPT::Pointer mpt, This *t): MPT::Listener(mpt), t_(t){}

    virtual void new_notification(typename MPT::Listener::Notification_type et){
      if (et == P::IS_EDITING && P::notifier()->is_editing()==false){
	if (P::notifier()->inserted_begin() != P::notifier()->inserted_end() 
	    || P::notifier()->erased_begin() != P::notifier()->erased_end()){
	  t_->update_tree();
	  t_->update_coordinates();
	}
      }
    }
  protected:
    This *t_;
  };
  friend class Table_listener;
  
private:
  //! This cannot be trivially copied with out ill effects
  Qt_moving_weighted_points_3(const This &){ CGAL_assertion(0);}
  const This operator=(const This &o) const {
    CGAL_assertion(0);
  }

public:
  
  //! This is for a bitfield
  typedef enum Draw_mode {POINT=1, SPHERE=0, TRANSPARENT_SPHERE=2} Draw_mode;

  //! Defaults to outline drawing
  Qt_moving_weighted_points_3( Traits tr,  
			       typename GUI::Pointer sim ): tr_(tr),
							    ik_(tr.instantaneous_kernel_object()), 
							    listener_(tr_.moving_point_table_pointer(), this),
							    guil_(sim, this),
							    siml_(sim->simulator(), this),
							    rt_(tr_.kinetic_kernel_object().reverse_time_object()){
   
    soss_= NULL;
    draw_labels_= true;
    point_size_= 5; 
    direction_of_time_=CGAL::POSITIVE;
    set_up_scene_graph(guil_.root());
  };

  virtual ~Qt_moving_weighted_points_3(){
    if (soss_!= NULL) {
      soss_->unschedule();
      delete soss_;
    }
  } 
  
  //! Access the SoCoordinate3 where it puts the current coordinates.
  SoCoordinate3 *coordinate_node() const {
    return coords_.get();
  }

  //void set_radius(double radius);
  void set_point_size(double ps);
  virtual void write(std::ostream &out) const;

  //! The draw mode.
  Draw_mode draw_mode(){
    return mode_;
  }

  //! Set whether points or balls are drawn
  void set_draw_mode(Draw_mode dm){
    mode_= dm;
    update_tree();
    update_coordinates();
  }

  //! Set whether the point labels are displayed
  bool draw_labels() const {
    return draw_labels_;
  }
  //! Set if the labels are drawn
  void set_draw_labels(bool tf) {
    draw_labels_=tf;
    update_tree();
    update_coordinates();
  }
protected:

  CGAL::Sign direction_of_time() const {
    return direction_of_time_;
  }
  void reverse_time();


  class Guil: public Gui_listener {
  public:
    Guil(typename GUI::Pointer& h, This *t): Gui_listener(h), t_(t){}
    void new_notification(typename Gui_listener::Notification_type nt){
      if (nt== Gui_listener::CURRENT_TIME){
	t_->update_coordinates();
      }
    }
  protected:
    This *t_;
  };
  friend class Guil;

  void set_up_scene_graph(SoSeparator* parent);

  void update_coordinates();

  void update_tree();

  double label_offset(){
    return .1;
  }
  double label_size(){
    return 2*.1;
  }

  static void update_callback(void *data, SoSensor *){
     This *th = reinterpret_cast<This*>(data);
     delete th->soss_;
     th->soss_=NULL;
     th->update_tree();
     th->update_coordinates();
  }

  

  static void keyboard_callback(void *data, SoEventCallback *eventCB){
    This *th = reinterpret_cast<This*>(data);
    const SoEvent *event= eventCB->getEvent();
    CGAL_assertion(event->isOfType(SoKeyboardEvent::getClassTypeId()));
    const SoKeyboardEvent *kbe= reinterpret_cast<const SoKeyboardEvent*>(event);
    //std::cout << "Pressed " << kbe->getPrintableCharacter() << std::endl;
    bool handled=false;
    if (kbe->getKey()== SoKeyboardEvent::L && kbe->getState()== SoButtonEvent::UP) {
      bool dl= th->draw_labels();
      handled= true;
      th->draw_labels_=!dl;
    } else if (kbe->getKey()== SoKeyboardEvent::S && kbe->getState()== SoButtonEvent::UP) {
      th->mode_=SPHERE;
      handled=true;
    } else if (kbe->getKey()== SoKeyboardEvent::P && kbe->getState()== SoButtonEvent::UP) {
      th->mode_=POINT;
      handled=true;
    } else if (kbe->getKey()== SoKeyboardEvent::T && kbe->getState()== SoButtonEvent::UP) {
      th->mode_=TRANSPARENT_SPHERE;
      handled=true;
    }
    if (handled){
      // let keystrokes be reused
      // eventCB->setHandled();
       if (th->soss_ ==NULL) {
	 th->soss_= new SoOneShotSensor(update_callback, th);
	 th->soss_->schedule();
       }
       
    }
  }
  
  unsigned int size() {
    /*unsigned int ct=0;
    for (typename MPT::Keys_iterator it= tr_.moving_point_table_pointer()->keys_begin();
	 it != tr_.moving_point_table_pointer()->keys_end(); ++it, ++ct);
	 return ct;*/
    return tr_.moving_point_table_pointer()->size();
  }

 
  Traits tr_;
  Draw_mode mode_;
  bool draw_labels_;
  double point_size_;
  IK ik_;
  //! I don't really want this mutable, but Inventor doesn't like constant nodes
  mutable Coin_pointer<SoCoordinate3> coords_;
  Coin_pointer<SoGroup> spheres_;
  Coin_pointer<SoShapeKit> points_;
  Coin_pointer<SoDrawStyle> style_;
  Coin_pointer<SoGroup> labels_;
  CGAL::Sign direction_of_time_;
  SoOneShotSensor* soss_;
  Table_listener listener_;
  Guil guil_;
  Siml siml_;
  typename Traits::Kinetic_kernel::Reverse_time rt_;
};


template <class T, class G>
void Qt_moving_weighted_points_3<T,G>::update_coordinates(){
  //std::cout << "updateing coordinates\n";
  //if (parent_==NULL) return;
  ik_.set_time(guil_.notifier()->current_time());
  
  coords_->point.setNum(size());
  SbVec3f *pts= coords_->point.startEditing();
  
  SbVec3f *vpts=NULL;
  if (points_!= NULL) {
    SoCoordinate3 *c= SO_GET_PART(points_, "coordinate3", SoCoordinate3);
    vpts= c->point.startEditing();
  }
  
  int cp=0;
  for (typename MPT::Keys_iterator it= tr_.moving_point_table_pointer()->keys_begin();
       it != tr_.moving_point_table_pointer()->keys_end(); ++it, ++cp){
    //std::cout << "drawing point " << *it  << "= " << ik_.to_static(*it) << std::endl;
    typename IK::Static_kernel::Weighted_point pt= ik_.static_object(*it);
    double w= CGAL::to_double(pt.weight());
    if (w < 0) w=0;
    double radius = std::sqrt(w);
    pts[it->index()].setValue(CGAL::to_double(pt.point().x()), CGAL::to_double(pt.point().y()),
			      CGAL::to_double(pt.point().z()));
    if (vpts != NULL) vpts[cp].setValue(CGAL::to_double(pt.point().x()),
					CGAL::to_double(pt.point().y()),
					CGAL::to_double(pt.point().z()));
    if (spheres_!= NULL){
      SoNode *n= spheres_->getChild(cp);
      CGAL_assertion(n->isOfType(SoShapeKit::getClassTypeId()));
      SoShapeKit *sh= reinterpret_cast<SoShapeKit*>(n);
      SoTransform *tr= SO_GET_PART(sh, "localTransform", SoTransform);
      tr->translation.setValue(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()), 
			       CGAL::to_double(pt.z()));
      SoSphere *sph= SO_GET_PART(sh, "shape", SoSphere);
      sph->radius.setValue(radius);
    } 
    if (labels_!= NULL){
      double offset;
      if (mode_ == POINT){
	offset=label_offset();
      } else {
	offset= 1.2*radius/std::sqrt(3.0);
      }
      SoNode *n= labels_->getChild(cp);
      CGAL_assertion(n->isOfType(SoShapeKit::getClassTypeId()));
      SoShapeKit *sh= reinterpret_cast<SoShapeKit*>(n);
      SoTransform *tr= SO_GET_PART(sh, "localTransform", SoTransform);
      tr->translation.setValue(CGAL::to_double(pt.point().x())+offset, CGAL::to_double(pt.point().y())+offset, 
			       CGAL::to_double(pt.point().z())+offset);
    } 
  }
  coords_->point.finishEditing();
  if (vpts!= NULL){
    SoCoordinate3 *c= SO_GET_PART(points_, "coordinate3", SoCoordinate3);
    c->point.finishEditing();
  }
}

  
template <class T, class G>
void Qt_moving_weighted_points_3<T,G>::update_tree(){
  //if (parent_==NULL) return; 
  int maxl=-1;
  int num=0;
  for (typename MPT::Keys_iterator it= tr_.moving_point_table_pointer()->keys_begin(); it != tr_.moving_point_table_pointer()->keys_end(); ++it){
    if (it->index() > maxl) maxl= it->index();
    ++num;
  }
  
  if (labels_ != NULL){
    guil_.root()->removeChild(labels_.get());
    labels_=NULL;
  }
  if (points_ != NULL){
    guil_.root()->removeChild(points_.get());
    points_=NULL;
  }
  if (spheres_ != NULL) {
    guil_.root()->removeChild(spheres_.get());
    spheres_=NULL;
  }
  if (maxl==-1) return;
  //std::cout << "updateing tree\n";
  coords_->point.setNum(maxl+1);
  if (mode_==POINT ){
    points_ = new SoShapeKit;
    Coin_pointer<SoCoordinate3> c= new SoCoordinate3;
    c->point.setNum(num);
    Coin_pointer<SoMaterial> mat= new SoMaterial;
    mat->diffuseColor.setValue(.8, 0,0);
    mat->ambientColor.setValue(.8,0,0);
    points_->setPart("material", mat.get());
    points_->setPart("coordinate3", c.get());
    Coin_pointer<SoPointSet> ps= new SoPointSet;
    ps->numPoints.setValue(num);
    points_->setPart("shape", ps.get());
    
    Coin_pointer<SoAppearanceKit> ak= new SoAppearanceKit;
    ak->setPart("drawStyle", style_.get());
    points_->setPart("appearance", ak.get());

    guil_.root()->addChild(points_.get());
  } else {
    spheres_= new SoGroup;
    guil_.root()->addChild(spheres_.get());
    Coin_pointer<SoMaterial> smat= new SoMaterial;
    smat->diffuseColor.setValue(.8, 0,0);
    if (mode_== TRANSPARENT_SPHERE){
      smat->transparency.setValue(.5);
    }
    for (int i=0; i< num; ++i){
      Coin_pointer<SoShapeKit> kit = new SoShapeKit;
      spheres_->addChild(kit.get());
      Coin_pointer<SoSphere> s= new SoSphere;
      s->radius.setValue(.01);
      kit->setPart("shape", s.get());
      Coin_pointer<SoTransform> tr= new SoTransform;
      kit->setPart("localTransform", tr.get());
      kit->setPart("material", smat.get());
    }
  }
  
  if (draw_labels_ != 0){
    Coin_pointer<SoMaterial> mat= new SoMaterial;
    mat->diffuseColor.setValue(1,1,1);
    mat->emissiveColor.setValue(1,1,1);
    labels_= new SoGroup;
    for (typename MPT::Keys_iterator kit = tr_.moving_point_table_pointer()->keys_begin();
	 kit != tr_.moving_point_table_pointer()->keys_end(); ++kit){
      Coin_pointer<SoShapeKit> k = new SoShapeKit;
      labels_->addChild(k.get());
      Coin_pointer<SoText2> s= new SoText2;
      std::string name = kit->string();
      s->string.setValue(name.c_str());
      k->setPart("shape", s.get());
      Coin_pointer<SoTransform> tr= new SoTransform;
      k->setPart("localTransform", tr.get());
      k->setPart("material", mat.get());
    }
    guil_.root()->addChild(labels_.get());
  } 
}

template <class T, class G>
void Qt_moving_weighted_points_3<T,G>::set_up_scene_graph(SoSeparator* parent){
  std::cout << "add to scene graph\n";
  SoEventCallback *myevcb= new SoEventCallback;
  myevcb->addEventCallback(SoKeyboardEvent::getClassTypeId(),keyboard_callback, this);
  parent->addChild(myevcb);

  style_=new SoDrawStyle;
  style_->pointSize.setValue(point_size_);
  
  coords_= new SoCoordinate3;
  
  parent->addChild(style_.get());
  parent->addChild(coords_.get());
  update_tree();
}



template <class T, class G>
void Qt_moving_weighted_points_3<T,G>::reverse_time(){
  //std::cout << "reversing time.\n";
  if (direction_of_time_== CGAL::POSITIVE) direction_of_time_=CGAL::NEGATIVE;
  else  direction_of_time_=CGAL::POSITIVE;

  tr_.moving_point_table_pointer()->set_is_editing(true);
  //typename MP::Traits::Reverse_time rt= tr_.moving_point_table_pointer()->traits_object().reverse_time_object();
  for (typename MPT::Keys_iterator kit= tr_.moving_point_table_pointer()->keys_begin(); kit != tr_.moving_point_table_pointer()->keys_end(); ++kit){
    tr_.moving_point_table_pointer()->set(*kit, rt_(tr_.moving_point_table_pointer()->at(*kit)));
  }
  tr_.moving_point_table_pointer()->set_is_editing(false);
}


//s->radius.setValue(radius_);

template <class T, class G>
void Qt_moving_weighted_points_3<T,G>::set_point_size(double ps){
  point_size_=ps;
  style_->pointSize.setValue(ps);
}

template <class T, class G>
void Qt_moving_weighted_points_3<T,G>::write(std::ostream &out) const {
  ik_.set_time(guil_.notifier()->current_time());
  for (typename MPT::Keys_iterator it= tr_.moving_point_table_pointer()->keys_begin(); 
       it != tr_.moving_point_table_pointer()->keys_end(); ++it){
    out << *it;
    out << ": " << ik_.static_object(*it) << std::endl;
  }
}
CGAL_KDS_END_NAMESPACE;

#endif // guard
