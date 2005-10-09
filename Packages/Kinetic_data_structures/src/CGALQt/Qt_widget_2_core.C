#include <CGAL/KDS/IO/internal/Qt_widget_2_core.h>

#include "Qt_widget_2_core.moc"

namespace CGAL {
  namespace KDS {
    namespace internal {
      void Qt_widget_2_core::redraw()  {
	  lock();
	  clear();
	  //std::cout << "size of drawables = " << drawable_s.size() << std::endl;
	  is_drawn_=false;
	  if (drawable_!= NULL) drawable_->new_notification(Listener::PICTURE_IS_CURRENT);
	  is_drawn_=true;
	  unlock();
	  //::CGAL::Qt_widget::redraw();
	}
      
      Qt_widget_2_core::Qt_widget_2_core(QMainWindow *parent): ::CGAL::Qt_widget(parent){
	drawable_=NULL;
	is_drawn_=false;
      }
    }
  }
}
