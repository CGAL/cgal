#ifndef CGAL_KDS_IO_INTERNAL_QT_TIMER_H
#define CGAL_KDS_IO_INTERNAL_QT_TIMER_H
#include <map>
#include <qtimer.h>
namespace CGAL{
  namespace KDS {
    namespace internal {
      class Qt_timer: QObject {
	Q_OBJECT
      public:
	class Listener{
	public:
	  typedef Qt_timer* Notifier_handle;
	  Listener(Notifier_handle h): h_(h){h->set_listener(this);}
	  typedef enum {TICKS} Notification_type;
	  virtual void new_notification(Notification_type) =0;
	  virtual ~Listener(){
	    h_->set_listener(NULL);
	  }
	protected:
	  Notifier_handle h_;
	};

	Qt_timer();

	int ticks() const {
	  return tick_;
	}
	void clear(){
	  //CGAL_precondition(id_!=-1);
	  if (id_!= -1) timer_.killTimer(id_);
	  id_=-1;
	};
	void run(double time_in_seconds);
      protected:
	QTimer timer_;
	Listener *cb_;
	int tick_;
	int id_;
	

	friend class Listener;
	void set_listener(Listener *l){
	  cb_=l;
	}

private slots:
void timerDone();
      };
    };
  };
};
#endif
