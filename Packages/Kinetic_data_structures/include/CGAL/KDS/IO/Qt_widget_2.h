#ifndef CGAL_KDS_QT_SIMULATOR_2_H_
#define CGAL_KDS_QT_SIMULATOR_2_H_
#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/IO/internal/GUI_base.h>
#include <CGAL/KDS/IO/internal/Qt_window_2.h>
#include <CGAL/KDS/IO/internal/Qt_widget_2_core.h>
#include <CGAL/KDS/IO/internal/Qt_timer.h>
#include <CGAL/KDS/Ref_counted.h>
#include <CGAL/KDS/Multi_listener.h>
#include <set>
#include <qapplication.h>

CGAL_KDS_BEGIN_NAMESPACE


//! A GUI for a kinetic simulation in 2D. 
/*!
  In order to draw graphics, extend Qt_gui_2::Listener. There is
  one notification type PICTURE_IS_VALID which is created
  whenever the validity of the picture changes, for example if
  the time has changed or the picture needs to be redrawn for
  some other reason. To redraw, use Qt_gui_2::Listener::widget()
  to get a CGAL::Qt_widget and use this widget to draw the
  desired geometry.

  As you might guess, Qt is required. 

  An example using this and Qt_moving_points_2 can be found in
  \example 2d_gui.cc.
*/
template <class Simulator_t >
class Qt_widget_2: public Ref_counted<Qt_widget_2<Simulator_t> > {
protected:
  typedef Qt_widget_2<Simulator_t> This;
  typedef Gui_base<Simulator_t, internal::Qt_timer> Graphical_base;
  typedef typename Simulator_t::Time Time;
  class Base_listener;
  typedef typename internal::Qt_core_listener<Graphical_base> Window_listener;
  typedef internal::Qt_widget_2_core Qt_widget;
  typedef internal::Qt_window_2 Qt_window;
  struct Listener_core{
    typedef enum {PICTURE_IS_VALID} Notification_type;
    typedef typename This::Pointer Notifier_pointer;
    Qt_widget* widget() const {return widget_;}
  private:
    friend class Qt_widget_2<Simulator_t>;
    void set_widget(Qt_widget *w){ widget_=w;}
  protected:
    Qt_widget *widget_;
  };
public:
  //! The base class for listeners for events.
  typedef Multi_listener<Listener_core> Listener;
  friend class  Multi_listener<Listener_core>;
  typedef Simulator_t Simulator;
  //typedef Const_ref_counted_pointer<This> Const_point;

  //! construct things
  Qt_widget_2(int argc, char *argv[],
	      typename  Simulator::Pointer sh,
	      double xmin=-10,double xmax=10, double ymin=-10, double ymax=10): app_(new QApplication(argc, argv)),
					     base_(new Graphical_base(sh)), 
					     base_l_( base_, this){
    
    app_->setMainWidget(new Qt_window(static_cast<int>(std::floor(xmin)),
				      static_cast<int>(std::ceil(xmax)), 
				      static_cast<int>(std::floor(ymin)), 
				      static_cast<int>(std::ceil(ymax))));
    window()->setCaption("KDS");
    window()->show();
    window_l_ = std::auto_ptr<Window_listener>(new Window_listener(window()->button_handler(), base_));
    widget_l_ = std::auto_ptr<Widget_listener>(new Widget_listener(widget(), this));
  }

  //! start the gui
  int begin_event_loop(){
    draw();
    return app_->exec();
  }
  
  //! Access a reference counted pointer to the simulator.
  /*!
    I am not sure that I need this method.
    \todo check if these methods are needed.
  */
  const typename Simulator::Pointer& simulator() const {
    bool let_me_know_if_this_is_used;
    return base_->simulator();
  }

  //! Access a reference counted pointer to the simulator.
  typename Simulator::Pointer& simulator() {
    return base_->simulator();
  }
  //! Return true if the current image of the scene is valid.
  /*!
    If this is false, then things need to redraw themselves.
  */
  bool picture_is_valid() const {
    return widget()->picture_is_current();
  }
  //! Return the current time as a double.
  double current_time() const {
    return base_->current_time();
  }
  //! Return the current speed of the simulation.
  double speed() const {
    return base_->speed();
  }
  //! Set the current speed of the simulation.
  void set_speed(double s) const {
    base_->set_speed(s);
  }
protected:

  //! Gui will call output_drawing
  void draw(){
    //CGAL_KDS_LOG(LOG_LOTS, "GUI: Drawing in gui.\n");
    for (typename std::set<Listener*>::iterator dit= drawables_.begin(); 
	 dit != drawables_.end(); ++dit){
      //log()->stream(Log::LOTS) << "GUI: Drawing something.\n";
      (*dit)->new_notification(Listener::PICTURE_IS_VALID);
    }
  }

  class Base_listener: public Graphical_base::Listener {
  public:
    Base_listener(typename Graphical_base::Pointer &b, This *t): Graphical_base::Listener(b), t_(t){}
    virtual void new_notification(typename Graphical_base::Listener::Notification_type nt){
      if (nt== Graphical_base::Listener::CURRENT_TIME){
	t_->widget()->set_picture_is_current(false);
      }
    }
  protected:
    This *t_;
  };

  class Widget_listener: public Qt_widget::Listener {
  public:
    Widget_listener(Qt_widget *w, This *t): Qt_widget::Listener(w),  t_(t){
    }
    virtual void new_notification(typename Qt_widget::Listener::Notification_type nt) {
      if (nt == Qt_widget::Listener::PICTURE_IS_CURRENT) {
	if (!widget()->picture_is_current()){
	  t_->draw();
	}
      }
    }
    virtual ~Widget_listener(){
    }
  protected:
    This *t_;
  };

  
  
private:
  void new_listener(Listener* d){
    CGAL_KDS_LOG(LOG_SOME, "GUI: Registered a drawable.\n");
    drawables_.insert(d);
    d->set_widget(widget());
    //widget()->set_picture_is_current(false);
  }
  void delete_listener(Listener* d){
    CGAL_KDS_LOG(LOG_SOME,"GUI: Unregistered a drawable.\n");
    drawables_.erase(d);
  }

  friend class Widget_listener; friend class Base_listener;

  Qt_window *window(){
    return reinterpret_cast<Qt_window*>(app_->mainWidget());
  }

  const Qt_window *window() const {
    return reinterpret_cast<const Qt_window*>(app_->mainWidget());
  }

  Qt_widget* widget() {
    return window()->widget();
  }
  
  const Qt_widget* widget() const {
    return window()->widget();
  }
  This operator=(const This &o) {
    return *this;
  }
  Qt_widget_2(const This &o){}
protected:
  std::auto_ptr<QApplication> app_;
  typename Graphical_base::Pointer base_;
  std::set<Listener*> drawables_;
  Base_listener base_l_;
  std::auto_ptr<Window_listener> window_l_;
  std::auto_ptr<Widget_listener> widget_l_;
};


CGAL_KDS_END_NAMESPACE

#endif // guard
