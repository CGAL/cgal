#include "Custom_manipulated_frame.h"

CustomManipulatedFrame::CustomManipulatedFrame(int timer_msec) : isHoldingMouse(false)
{
  timer = new QTimer(this);
  connect(timer, SIGNAL(timeout()), this, SLOT(timerUpdate()));
  timer->start(timer_msec);
}

void CustomManipulatedFrame::mousePressEvent(QMouseEvent* const event, qglviewer::Camera* const camera)
{
  isHoldingMouse = true;
  ManipulatedFrame::mousePressEvent(event, camera);
}

void CustomManipulatedFrame::mouseReleaseEvent(QMouseEvent* const event, qglviewer::Camera* const camera)
{
  isHoldingMouse = false;
  ManipulatedFrame::mouseReleaseEvent(event, camera);
}

void CustomManipulatedFrame::timerUpdate()
{
  if(isHoldingMouse)
  { emit modified(); }
}

#include "Custom_manipulated_frame.moc"
