// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef CAMERA_MANIP_ROT_H
#define CAMERA_MANIP_ROT_H

#include <qevent.h>
#include <qvector2d.h>

#include "Camera_manip.h"

class Camera_manip_rot : public Camera_manip {
public:
  Camera_manip_rot(Camera& camera);

protected:
  virtual void mouse_move_event(QMouseEvent* e) override;

private:
  float m_theta = 0;
  float m_phi = 0;
};

#endif
