// Copyright(c) 2012, 2020  Tel - Aviv University(Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef SINGLE_VERTEX_H
#define SINGLE_VERTEX_H

#include <vector>
#include <qvector3d.h>
#include "Common_defs.h"


class SingleVertex : protected OpenGLFunctionsBase
{
public:
  SingleVertex(const QVector3D& pos);

  void set_visible(bool flag);
  void set_pos(const QVector3D& pos);
  const QVector3D& get_pos() const;

  void draw();

private:
  bool      m_visible;
  GLuint    m_vao, m_vbo;
  QVector3D m_pos;
};


#endif
