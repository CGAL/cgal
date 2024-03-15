// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef CAMERA_MANIP_H
#define CAMERA_MANIP_H

#include <qevent.h>
#include <qvector2d.h>

#include "Camera.h"
#include "GUI_event_handler.h"

class Camera_manip : public GUI_event_handler {
public:
  Camera_manip(Camera& camera);

protected:
  Camera& m_camera;
};


#endif
