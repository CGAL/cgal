#include <CGAL/KDS/IO/internal/Qt_timer.h>

#include "Qt_timer.moc"

namespace CGAL {
  namespace KDS {
    namespace internal {
      
      Qt_timer::Qt_timer(): tick_(0), id_(-1){
	connect( &timer_, SIGNAL(timeout()),
		 this, SLOT(timerDone()) );
	
      }

      void Qt_timer::run(double time_in_seconds){
	id_=timer_.start(static_cast<int>(time_in_seconds*1000), true);
	//CGAL_postcondition(id_ != -1);
      }


      void Qt_timer::timerDone(){
	++tick_;
	cb_->new_notification(Listener::TICKS);
      }
    }
  }
}
