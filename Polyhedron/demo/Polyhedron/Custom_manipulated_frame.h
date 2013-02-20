#ifndef CUSTOM_MANIPULATED_FRAME_H
#define CUSTOM_MANIPULATED_FRAME_H
#include "Custom_manipulated_frame_config.h"

#include <QObject>
#include <QTimer>
#include <QGLViewer/manipulatedFrame.h>
/*
 * Holds mouse pressing state by interleaving mousePressEvent and mouseReleaseEvent of base class --ManipulatedFrame.
 * Uses a timer to emit modified signal of base class while mouse is pressed.
 */
class CUSTOM_MANIPULATED_FRAME_EXPORT CustomManipulatedFrame 
  : public qglviewer::ManipulatedFrame 
{
	Q_OBJECT
	public:
	  CustomManipulatedFrame(int timer_msec = 50/*20 fps*/);
	protected:
	  virtual void mousePressEvent(QMouseEvent* const event, qglviewer::Camera* const camera);
	  virtual void mouseReleaseEvent(QMouseEvent* const event, qglviewer::Camera* const camera);

	public slots:
	  void timerUpdate();

	private:
	  bool isHoldingMouse; 
	  QTimer* timer;
};
#endif