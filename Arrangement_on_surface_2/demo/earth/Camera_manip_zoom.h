
#ifndef CAMERA_MANIP_ZOOM_H
#define CAMERA_MANIP_ZOOM_H

#include <qevent.h>
#include <qvector2d.h>

#include "Camera_manip.h"


class Camera_manip_zoom : public Camera_manip
{
public:
  Camera_manip_zoom(Camera& camera);

protected:
  virtual void mouse_move_event(QMouseEvent* e) override;

private:
};

#endif
