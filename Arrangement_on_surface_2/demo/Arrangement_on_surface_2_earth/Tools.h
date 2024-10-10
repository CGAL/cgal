// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef TOOLS_H
#define TOOLS_H

#include <string>

#include <QVector2D>
#include <QVector3D>

std::string read_file(const std::string& file_name);

std::ostream& operator << (std::ostream& os, const QVector2D& v);
std::ostream& operator << (std::ostream& os, const QVector3D& v);
std::ostream& operator << (std::ostream& os, const QVector4D& v);


#endif
