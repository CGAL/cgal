#ifndef CGAL_KDS_IO_MOVING_POINT_SET_2_H
#define CGAL_KDS_IO_MOVING_POINT_SET_2_H
#include <CGAL/KDS/basic.h>
//#include <CGAL/KDS/Moving_object_table.h>
#include <CGAL/KDS/Cartesian_instantaneous_kernel.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/KDS/Ref_counted.h>
#include <CGAL/KDS/Simulator_objects_listener.h>


CGAL_KDS_BEGIN_NAMESPACE


//! A graphical moving point set
/*!  This allows a set of moving points to be used in a
  CGAL::KDS::Qt_gui_2 GUI. It listens for PICTURE_IS_VALID
  notifications from the CGAL::KDS::Qt_gui_2 and DIRECTION_OF_TIME
  notifcations from the CGAL::KDS::Simulator and handles them
  accordingly.
 
*/
template <class TraitsC, class GUI> 
class Qt_moving_points_2: public Ref_counted<Qt_moving_points_2<TraitsC, GUI> >  {
protected:
  typedef TraitsC Traits;
  //typedef Moving_object_set<Moving_point_t> Parent;
  typedef Qt_moving_points_2<Traits, GUI> This;

  typedef typename Traits::Simulator Simulator;
  typedef typename GUI::Listener Gui_listener;
  typedef typename Simulator::Listener Simulator_listener;
  typedef CGAL::KDS::Simulator_objects_listener<Simulator_listener, This> Siml;
  friend class CGAL::KDS::Simulator_objects_listener<Simulator_listener, This>;
  class Guil;
public:
  //typedef typename CGAL::Ref_counted_pointer<This> Pointer;
  //typedef typename CGAL::Const_ref_counted_pointer<This> Const_handle;

  //! Only outline drawing is currently supported
  typedef enum Draw_mode {FILLED, POINT, OUTLINE} Draw_mode;

  //! Defaults to outline drawing
  Qt_moving_points_2(typename GUI::Pointer sim, 
		     Traits tr): traits_(tr),
				 ik_(tr.instantaneous_kernel_object()),
				 _mode(POINT), _radius(.2),
				 direction_of_time_(CGAL::POSITIVE),
				 guil_(sim, const_cast<This*>(this)), 
				 siml_(tr.simulator_pointer(),const_cast<This*>( this)),
				 rt_(tr.kinetic_kernel_object().reverse_time_object())
  {
  };

  virtual ~Qt_moving_points_2(){}

  //! Change how things are drawn
  void set_draw_mode(Draw_mode mode){
    _mode=mode;
  }

  //! change how big things are drawn
  void set_radius(double radius){
    _radius=radius;
  }
  
  //! Set the CGAL::Sign field corresponding to the direction of time.
  /*!
    This can have the effect of reversing the motions of the points. 
  */
  CGAL::Sign direction_of_time() const {
    return direction_of_time_;
  }

  /*virtual void write(std::ostream &out) const {
    ik_.set_time(sim_->current_time_nt());
    for (typename Moving_point_table::key_iterator it= Moving_point_table::begin_keys(); it != Moving_point_table::end_keys(); ++it){
    out << *it;
    out << ": " << ik_.to_static(*it) << std::endl;
    }
    }*/

protected:
  
  class Guil: public Gui_listener {
    typedef Gui_listener P;
  public:
    Guil(typename GUI::Pointer &h, This *t): Gui_listener(h), t_(t){}
    void new_notification(typename Gui_listener::Notification_type nt){
      if (nt== Gui_listener::PICTURE_IS_VALID && !P::notifier()->picture_is_valid()){
	t_->draw();
      }
    }
  protected:
    This *t_;
  };
  friend class Guil;

  void draw() const;

  void reverse_time();

  Traits traits_;
  typename Traits::Instantaneous_kernel ik_;
  Draw_mode _mode;
  double _radius;
  CGAL::Sign direction_of_time_;
  Guil guil_;
  Siml siml_;
  typename Traits::Kinetic_kernel::Reverse_time rt_;
};

template <class T, class G>
void Qt_moving_points_2<T,G>::draw() const {
  //std::cout << "Drawing mpt MPT\n";
  typedef typename Traits::Static_kernel::Point_2 P2;
  typedef typename Traits::Static_kernel::Circle_2 C;
  typename Traits::NT ntt(guil_.notifier()->current_time());
  ik_.set_time(ntt);
  
  CGAL::Qt_widget *w= guil_.widget();
  *w << CGAL::PointSize(10) << CGAL::LineWidth(1);
  *w << CGAL::FillColor(CGAL::Color(0,0,0));
  //out.setFillColor(CGAL::Color(0,255,0));
  *w << CGAL::Color(0,0,0);
  //out << C(P2(0,0), 2) << C(P2(0,0), 1);
  //out << CGAL::BackgroundColor(CGAL::Color(125,125,125));
  for (typename Traits::Moving_point_table::Keys_iterator it= traits_.moving_point_table_pointer()->keys_begin();
       it != traits_.moving_point_table_pointer()->keys_end(); ++it){
    //std::cout << "drawing point " << *it  << "= " << ik_.to_static(*it) << std::endl;
    if (_mode== OUTLINE) {
      *w << C(ik_.static_object(*it), _radius);
    } else if (_mode == POINT){
      *w << ik_.static_object(*it);
    }
  }
}

template <class T, class G>
void Qt_moving_points_2<T,G>::reverse_time(){
  if (direction_of_time_ == CGAL::POSITIVE) direction_of_time_=CGAL::NEGATIVE;
  else direction_of_time_= CGAL::POSITIVE;
  
  traits_.moving_point_table_pointer()->set_is_editing(true);
  for (typename Traits::Moving_point_table::Keys_iterator kit= traits_.moving_point_table_pointer()->keys_begin(); 
       kit != traits_.moving_point_table_pointer()->keys_end(); ++kit){
    traits_.moving_point_table_pointer()->set(*kit, rt_(traits_.moving_point_table_pointer()->at(*kit)));
  }
  traits_.moving_point_table_pointer()->set_is_editing(false);
}


CGAL_KDS_END_NAMESPACE

#endif // guard
