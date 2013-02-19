#include "Custom_manipulated_frame.h"

CustomManipulatedFrame::CustomManipulatedFrame() : isHoldingMouse(false)
{
  std::cerr << "-----------------------------------------" << std::endl;
  std::cerr << "Construct" << std::endl;
  std::cerr << "-----------------------------------------" << std::endl;
  timer = new QTimer(this);
  connect(timer, SIGNAL(timeout()), this, SLOT(timerUpdate()));
  timer->start(1000);
}

void CustomManipulatedFrame::mousePressEvent(QMouseEvent* const event, qglviewer::Camera* const camera)
{
  std::cerr << "-----------------------------------------" << std::endl;
  std::cerr << "Holding TRUE" << std::endl;
  std::cerr << "-----------------------------------------" << std::endl;
  isHoldingMouse = true;
  ManipulatedFrame::mousePressEvent(event, camera);
}

void CustomManipulatedFrame::mouseReleaseEvent(QMouseEvent* const event, qglviewer::Camera* const camera)
{
  std::cerr << "-----------------------------------------" << std::endl;
  std::cerr << "Holding FALSE" << std::endl;
  std::cerr << "-----------------------------------------" << std::endl;
  isHoldingMouse = false;
  ManipulatedFrame::mouseReleaseEvent(event, camera);
}

void CustomManipulatedFrame::timerUpdate()
{
  std::cerr << "Timer UPDATE "<<isHoldingMouse << std::endl;
  if(isHoldingMouse)
    Q_EMIT manipulated();
}

#include "Custom_manipulated_frame.moc"
