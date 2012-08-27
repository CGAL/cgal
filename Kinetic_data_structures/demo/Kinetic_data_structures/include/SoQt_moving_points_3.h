// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_QT_MOVING_POINT_TABLE_3_H
#define CGAL_KINETIC_QT_MOVING_POINT_TABLE_3_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/Simulator_objects_listener.h>
#include "SoQt_handle.h"

#include <Inventor/events/SoButtonEvent.h>
#include <Inventor/events/SoEvent.h>
#include <Inventor/events/SoKeyboardEvent.h>
#include <Inventor/nodekits/SoAppearanceKit.h>
#include <Inventor/nodekits/SoShapeKit.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoPointSet.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoText2.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/sensors/SoOneShotSensor.h>

namespace CGAL { namespace Kinetic {;

//! A graphical moving point set in 3D
/*!
  This class shows how to interact with the Coin and Qt based 3D gui.

  The appearance can be controlled from the keyboard.
  - L toggles labels
  - S draws the points as constant radius spheres
  - P draws the points as points
*/
template <class SimTraits, class GUI>
class SoQt_moving_points_3: public Ref_counted< SoQt_moving_points_3<SimTraits, GUI> >
{
protected:
  typedef SoQt_moving_points_3<SimTraits, GUI> This;
  typedef typename SimTraits::Active_points_3_table MPT;
  typedef typename SimTraits::Simulator Simulator;
  typedef typename GUI::Listener Gui_listener;
  typedef typename SimTraits::Instantaneous_kernel IK;
  typedef typename MPT::Data::Coordinate::NT NT;

 
  typedef typename SimTraits::Active_points_3_table Moving_point_table;

  CGAL_KINETIC_LISTEN1(Simulator, DIRECTION_OF_TIME, reverse_time());
  CGAL_KINETIC_LISTEN1(GUI, CURRENT_TIME, update_coordinates());
  CGAL_KINETIC_LISTEN1(MPT, IS_EDITING, editing());
  void editing() {
    if (CGAL_KINETIC_NOTIFIER(MPT)->inserted_begin() != CGAL_KINETIC_NOTIFIER(MPT)->inserted_end()
        || CGAL_KINETIC_NOTIFIER(MPT)->erased_begin() != CGAL_KINETIC_NOTIFIER(MPT)->erased_end()) {
        update_tree();
        update_coordinates();
      }
  }
private:
  //! This cannot be trivially copied with out ill effects
  SoQt_moving_points_3(const This &){CGAL_error();}
  const This operator=(const This &o) const
  {

    CGAL_error();
  }

public:
  //! The different modes for drawing
  typedef enum Draw_mode {POINT=1, SPHERE=0}
    Draw_mode;

  //! Defaults to outline drawing
  SoQt_moving_points_3(SimTraits tr, typename GUI::Handle sim):  tr_(tr),
                                                                 ik_(tr.instantaneous_kernel_object()),
                                                                 rt_(tr.kinetic_kernel_object().reverse_time_object()) {
    soss_= NULL;
    draw_labels_= true;
    radius_=.1;
    point_size_= 5;
    direction_of_time_=CGAL::POSITIVE;
    sim->soqt_examiner_viewer_pointer()->getPointSizeLimits(point_size_bounds_, point_size_granularity_);
    CGAL_KINETIC_INIT_LISTEN(MPT, tr_.active_points_3_table_handle());
    CGAL_KINETIC_INIT_LISTEN(Simulator, tr.simulator_handle());
    CGAL_KINETIC_INIT_LISTEN(GUI, sim);
    set_up_scene_graph(root());
  };

  virtual ~SoQt_moving_points_3() {
    if (soss_!= NULL) {
      soss_->unschedule();
      delete soss_;
    }
  }

  //! The node where the coordinates are put.
  SoCoordinate3 *coordinate_node() const
  {
    return coords_.get();
  }

  //! Set the radius to draw the points.
  void set_radius(double radius);
  //! The the point size.
  void set_point_size(double ps);
  virtual void write(std::ostream &out) const;

  //! The current draw mode
  Draw_mode draw_mode() {
    return mode_;
  }

  //! Set whether points or balls are drawn
  void set_draw_mode(Draw_mode dm) {
    mode_= dm;
    update_tree();
    update_coordinates();
  }

  //! Set whether the point labels are displayed
  bool draw_labels() const
  {
    return draw_labels_;
  }
  //! Set the draw_labels bit.
  void set_draw_labels(bool tf) {
    draw_labels_=tf;
    update_tree();
    update_coordinates();
  }
protected:
  SoSeparator* root() const {
    return listener_GUI_.root();
  }
  CGAL::Sign direction_of_time() const
  {
    return direction_of_time_;
  }
  void reverse_time();

  class Guil: public Gui_listener
  {
  public:
    Guil(typename GUI::Handle& h, This *t): Gui_listener(h), t_(t){}
    void new_notification(typename Gui_listener::Notification_type nt) {
      if (nt== Gui_listener::CURRENT_TIME) {
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

  double label_offset() {
    return radius_;
  }
  double label_size() {
    return 2*radius_;
  }

  static void update_callback(void *data, SoSensor *) {
    This *th = reinterpret_cast<This*>(data);
    delete th->soss_;
    th->soss_=NULL;
    th->update_tree();
    th->update_coordinates();
  }

  static void keyboard_callback(void *data, SoEventCallback *eventCB) {
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
    }
    else if (kbe->getKey()== SoKeyboardEvent::S && kbe->getState()== SoButtonEvent::UP) {
      th->mode_=SPHERE;
      handled=true;
    }
    else if (kbe->getKey()== SoKeyboardEvent::P && kbe->getState()== SoButtonEvent::UP) {
      th->mode_=POINT;
      handled=true;
    }
    if (handled) {
      //eventCB->setHandled();
      if (th->soss_ ==NULL) {
	th->soss_= new SoOneShotSensor(update_callback, th);
	th->soss_->schedule();
      }
    }
  }

  unsigned int size() {
    unsigned int ct=0;
    for (typename MPT::Key_iterator it= tr_.active_points_3_table_handle()->keys_begin();
	 it != tr_.active_points_3_table_handle()->keys_end(); ++it, ++ct);
    return ct;
  }

  float fix_point_size(double ps) const {
    if (ps <= point_size_bounds_[0]) return point_size_bounds_[0];
    else if (ps >= point_size_bounds_[1]) return point_size_bounds_[1];
    return ps;
  }

  Draw_mode mode_;
  bool draw_labels_;
  double radius_;
  SbVec2f point_size_bounds_;
  float point_size_granularity_;
  double point_size_;
  SimTraits tr_;
  IK ik_;
  //! I don't really want this mutable, but Inventor doesn't like constant nodes
  mutable SoQt_handle<SoCoordinate3> coords_;
  SoQt_handle<SoGroup> spheres_;
  SoQt_handle<SoShapeKit> points_;
  SoQt_handle<SoDrawStyle> style_;
  SoQt_handle<SoGroup> labels_;
  CGAL::Sign direction_of_time_;
  SoOneShotSensor* soss_;
  typename SimTraits::Kinetic_kernel::Reverse_time rt_;
};

template <class Tr, class G>
void SoQt_moving_points_3<Tr, G>::update_coordinates()
{
  //std::cout << "updateing coordinates\n";
  //if (parent_==NULL) return;
  ik_.set_time(NT(CGAL_KINETIC_NOTIFIER(GUI)->current_time()));

  coords_->point.setNum(size());
  SbVec3f *pts= coords_->point.startEditing();

  SbVec3f *vpts=NULL;
  if (points_!= NULL) {
    SoCoordinate3 *c= SO_GET_PART(points_, "coordinate3", SoCoordinate3);
    vpts= c->point.startEditing();
  }

  int cp=0;
  for (typename MPT::Key_iterator it= tr_.active_points_3_table_handle()->keys_begin();
       it != tr_.active_points_3_table_handle()->keys_end(); ++it, ++cp) {
    //std::cout << "drawing point " << *it  << "= " << ik_.to_static(*it) << std::endl;
    typename IK::Static_kernel::Point_3 pt= ik_.current_coordinates_object()(*it);
    pts[it->index()].setValue(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()),
			      CGAL::to_double(pt.z()));
    if (vpts != NULL) vpts[cp].setValue(CGAL::to_double(pt.x()),
					CGAL::to_double(pt.y()),
					CGAL::to_double(pt.z()));
    if (spheres_!= NULL) {
      SoNode *n= spheres_->getChild(cp);
      CGAL_assertion(n->isOfType(SoShapeKit::getClassTypeId()));
      SoShapeKit *sh= reinterpret_cast<SoShapeKit*>(n);
      SoTransform *tr= SO_GET_PART(sh, "localTransform", SoTransform);
      tr->translation.setValue(CGAL::to_double(pt.x()), CGAL::to_double(pt.y()),
			       CGAL::to_double(pt.z()));
    }
    if (labels_!= NULL) {
      SoNode *n= labels_->getChild(cp);
      CGAL_assertion(n->isOfType(SoShapeKit::getClassTypeId()));
      SoShapeKit *sh= reinterpret_cast<SoShapeKit*>(n);
      SoTransform *tr= SO_GET_PART(sh, "localTransform", SoTransform);
      tr->translation.setValue(CGAL::to_double(pt.x())+radius_, CGAL::to_double(pt.y())+radius_,
			       CGAL::to_double(pt.z())+radius_);
    }
  }
  coords_->point.finishEditing();
  if (vpts!= NULL) {
    SoCoordinate3 *c= SO_GET_PART(points_, "coordinate3", SoCoordinate3);
    c->point.finishEditing();
  }
}


template <class Tr, class G>
void SoQt_moving_points_3<Tr, G>::update_tree()
{
  //if (parent_==NULL) return;
  int maxl=-1;
  int num=0;
  for (typename MPT::Key_iterator it= tr_.active_points_3_table_handle()->keys_begin(); it != tr_.active_points_3_table_handle()->keys_end(); ++it) {
    if (static_cast<int>(it->index()) > maxl) maxl= it->index();
    ++num;
  }

  if (labels_ != NULL) {
    root()->removeChild(labels_.get());
    labels_=NULL;
  }
  if (points_ != NULL) {
    root()->removeChild(points_.get());
    points_=NULL;
  }
  if (spheres_ != NULL) {
    root()->removeChild(spheres_.get());
    spheres_=NULL;
  }
  if (maxl==-1) return;
  //std::cout << "updateing tree\n";
  coords_->point.setNum(maxl+1);
  if (mode_==POINT ) {
    points_ = new SoShapeKit;
    SoQt_handle<SoCoordinate3> c= new SoCoordinate3;
    c->point.setNum(num);
    SoQt_handle<SoMaterial> mat= new SoMaterial;
    mat->diffuseColor.setValue(.8, 0,0);
    mat->ambientColor.setValue(.8,0,0);
    points_->setPart("material", mat.get());
    points_->setPart("coordinate3", c.get());
    SoQt_handle<SoPointSet> ps= new SoPointSet;
    ps->numPoints.setValue(num);
    points_->setPart("shape", ps.get());

    SoQt_handle<SoAppearanceKit> ak= new SoAppearanceKit;
    ak->setPart("drawStyle", style_.get());
    points_->setPart("appearance", ak.get());

    root()->addChild(points_.get());
  }
  else {
    spheres_= new SoGroup;
    root()->addChild(spheres_.get());
    SoQt_handle<SoMaterial> smat= new SoMaterial;
    smat->diffuseColor.setValue(.8, 0,0);
    for (int i=0; i< num; ++i) {
      SoQt_handle<SoShapeKit> kit = new SoShapeKit;
      spheres_->addChild(kit.get());
      SoQt_handle<SoSphere> s= new SoSphere;
      s->radius.setValue(radius_);
      kit->setPart("shape", s.get());
      SoQt_handle<SoTransform> tr= new SoTransform;
      kit->setPart("localTransform", tr.get());
      kit->setPart("material", smat.get());
    }
  }

  if (draw_labels_ != 0) {
    SoQt_handle<SoMaterial> mat= new SoMaterial;
    mat->diffuseColor.setValue(1,1,1);
    mat->emissiveColor.setValue(1,1,1);
    labels_= new SoGroup;
    for (typename MPT::Key_iterator kit = tr_.active_points_3_table_handle()->keys_begin();
	 kit != tr_.active_points_3_table_handle()->keys_end(); ++kit) {
      SoQt_handle<SoShapeKit> k = new SoShapeKit;
      labels_->addChild(k.get());
      SoQt_handle<SoText2> s= new SoText2;
      std::string name = kit->string();
      s->string.setValue(name.c_str());
      k->setPart("shape", s.get());
      SoQt_handle<SoTransform> tr= new SoTransform;
      k->setPart("localTransform", tr.get());
      k->setPart("material", mat.get());
    }
    root()->addChild(labels_.get());
  }
}


template <class Tr, class G>
void SoQt_moving_points_3<Tr, G>::set_up_scene_graph(SoSeparator* parent)
{
  std::cout << "add to scene graph\n";
  SoEventCallback *myevcb= new SoEventCallback;
  myevcb->addEventCallback(SoKeyboardEvent::getClassTypeId(),keyboard_callback, this);
  parent->addChild(myevcb);

  style_=new SoDrawStyle;
  style_->pointSize.setValue(fix_point_size(point_size_));

  coords_= new SoCoordinate3;

  parent->addChild(style_.get());
  parent->addChild(coords_.get());
  update_tree();
}


template <class Tr, class G>
void SoQt_moving_points_3<Tr, G>::reverse_time()
{
  //std::cout << "reversing time.\n";
  if (direction_of_time_== CGAL::POSITIVE) direction_of_time_=CGAL::NEGATIVE;
  else  direction_of_time_=CGAL::POSITIVE;

  tr_.active_points_3_table_handle()->set_is_editing(true);
  //typename MP::Traits::Reverse_time rt= tr_.active_points_3_table_pointer()->traits_object().reverse_time_object();
  for (typename Tr::Active_points_3_table::Key_iterator kit= tr_.active_points_3_table_handle()->keys_begin(); kit != tr_.active_points_3_table_handle()->keys_end(); ++kit) {
    tr_.active_points_3_table_handle()->set(*kit, rt_(tr_.active_points_3_table_handle()->at(*kit)));
  }
  tr_.active_points_3_table_handle()->set_is_editing(false);
}


template <class Tr, class G>
void SoQt_moving_points_3<Tr, G>:: set_radius(double radius)
{
  radius_=radius;
  if (spheres_!= NULL) {
    for (int i=0; i< spheres_->getNumChildren(); ++i) {
      SoNode *c= spheres_->getChild(i);
      CGAL_assertion(c->isOfType(SoShapeKit::getClassTypeId()));
      SoShapeKit *kit= reinterpret_cast<SoShapeKit*>(c);
      SoSphere *s= SO_GET_PART(kit, "shape", SoSphere);
      s->radius.setValue(radius_);
    }
  }
}


template <class Tr, class G>
void SoQt_moving_points_3<Tr, G>::set_point_size(double ps)
{
  point_size_=ps;
  style_->pointSize.setValue(ps);
}


template <class Tr, class G>
void SoQt_moving_points_3<Tr, G>::write(std::ostream &out) const
{
  ik_.set_time(NT(CGAL_KINETIC_NOTIFIER(GUI)->current_time()));
  for (typename MPT::Key_iterator it= tr_.active_points_3_table_handle()->keys_begin();
       it != tr_.active_points_3_table_handle()->keys_end(); ++it) {
    out << *it;
    out << ": " << ik_.current_coordinates_object()(*it) << std::endl;
  }
}


} } //namespace CGAL::Kinetic;
#endif                                            // guard
