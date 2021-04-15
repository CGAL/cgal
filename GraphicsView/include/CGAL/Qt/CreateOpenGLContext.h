// Copyright (c) 2015  GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau and Maxime Gimeno
#ifndef CGAL_QT_CREATE_OPENGL_CONTEXT_H
#define CGAL_QT_CREATE_OPENGL_CONTEXT_H

#include <QOpenGLContext>
#include <QGLContext>
namespace CGAL{
namespace Qt{
inline QGLContext* createOpenGLContext()
{
    QOpenGLContext *context = new QOpenGLContext();
    QSurfaceFormat format;
    format.setVersion(2,1);
    format.setProfile(QSurfaceFormat::CompatibilityProfile);
    context->setFormat(format);
    QGLContext *result = QGLContext::fromOpenGLContext(context);
    result->create();
    return result;
}
} // namespace Qt
} // namespace CGAL
#endif
