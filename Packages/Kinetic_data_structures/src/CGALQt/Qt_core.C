#include <CGAL/KDS/IO/internal/Qt_core.h>

#include "Qt_core.moc"

namespace CGAL {
  namespace KDS {
    namespace internal {
      Qt_core::Qt_core(){
	playable_=NULL;
	  mode_=STOP;
      }
   
      void Qt_core::play_button(){
	CGAL_KDS_LOG(LOG_SOME, "Play button pushed.\n");
	mode_=RUN;
	if (playable_ != NULL) {
	  playable_->new_notification(Listener::LAST_BUTTON_PRESSED);
	} else {
	  CGAL_KDS_ERROR("...but no handler was registered.\n");
	}
      };

      void Qt_core::pause_button(){
	CGAL_KDS_LOG(LOG_SOME, "Pause button pushed.\n");
	mode_=PAUSE;
	if (playable_ != NULL) {
	  playable_->new_notification(Listener::LAST_BUTTON_PRESSED);
	} else {
	  CGAL_KDS_ERROR("...but no handler was registered.\n");
	}
      }

      void Qt_core::stop_button(){
	CGAL_KDS_LOG(LOG_SOME, "Stop button pushed.\n");
	mode_=STOP;
	if (playable_ != NULL) {
	  playable_->new_notification(Listener::LAST_BUTTON_PRESSED);
	} else {
	  CGAL_KDS_ERROR( "...but no handler was registered.\n");
	}
      }

      void Qt_core::play_to_button(){
	CGAL_KDS_LOG(LOG_SOME, "Play_to button pushed.\n");
	mode_=RUN_TO;
	if (playable_ != NULL) {
	  playable_->new_notification(Listener::LAST_BUTTON_PRESSED);
	} else {
	  CGAL_KDS_ERROR( "...but no handler was registered.\n");
	}
      }

      void Qt_core::play_through_button(){
	CGAL_KDS_LOG(LOG_SOME, "Play through button pushed.\n");
	mode_= RUN_THROUGH;
	if (playable_ != NULL) {
	  playable_->new_notification(Listener::LAST_BUTTON_PRESSED);
	} else {
	  CGAL_KDS_ERROR( "...but no handler was registered.\n");
	}
      }
      void Qt_core::reverse_button(){
	CGAL_KDS_LOG(LOG_SOME, "Reverse button pushed.\n");
	mode_=REVERSE;
	if (playable_ != NULL) {
	  playable_->new_notification(Listener::LAST_BUTTON_PRESSED);
	} else {
	  CGAL_KDS_ERROR("...but no handler was registered.\n");
	}
      }
      void Qt_core::faster_button(){
	CGAL_KDS_LOG(LOG_SOME, "Faster button pushed.\n");
	mode_=FASTER;
	if (playable_ != NULL) {
	  playable_->new_notification(Listener::LAST_BUTTON_PRESSED);
	} else {
	  CGAL_KDS_ERROR( "...but no handler was registered.\n");
	}
      }
      void Qt_core::slower_button(){
	CGAL_KDS_LOG(LOG_SOME, "Slower button pushed.\n");
	mode_=SLOWER;
	if (playable_ != NULL) {
	  playable_->new_notification(Listener::LAST_BUTTON_PRESSED);
	} else {
	  CGAL_KDS_ERROR("...but no handler was registered.\n");
	}
      }
    }
  }
}
