// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef WORLD_COORD_AXES_H
#define WORLD_COORD_AXES_H

#include <qvector4d.h>

#include "Common_defs.h"

class World_coord_axes : protected OpenGLFunctionsBase {
public:
  World_coord_axes(float length);

  void draw();

private:
  GLuint m_vao;
  GLuint m_vbo;
};


#endif
