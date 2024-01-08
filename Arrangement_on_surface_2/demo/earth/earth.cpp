// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#include <boost/filesystem.hpp>

#include <QApplication>
#include <QLabel>
#include <QSurfaceFormat>

#ifndef QT_NO_OPENGL
#include "Main_widget.h"
#endif

namespace fi = boost::filesystem;

int main(int argc, char* argv[]) {
  QApplication app(argc, argv);
  auto args = QCoreApplication::arguments();
  if (args.size() > 2) {
    qDebug() << "Usage: earth [<arragement-file.json>]";
    return(-1);
  }
  const char* file_name = (argc > 1) ? argv[1] :
    "../../../Data/data/arrangements_3/sphere/ne_110m_admin_0_countries.json";
  fi::path file_path(file_name);
  if (! fi::exists(file_path)) {
    std::cerr << "Error: failed to find file " << file_name << "\n";
    return -1;
  }

  QSurfaceFormat format;
  format.setVersion(3, 3);
  format.setProfile(QSurfaceFormat::CoreProfile);
  //format.setProfile(QSurfaceFormat::CompatibilityProfile);
  //format.setOptions(QSurfaceFormat::DeprecatedFunctions);
  //QSurfaceFormat::setDefaultFormat(format);
  format.setDepthBufferSize(24);
  QSurfaceFormat::setDefaultFormat(format);

  app.setApplicationName("Earth");
  app.setApplicationVersion("0.1");
  try {
#ifndef QT_NO_OPENGL
    Main_widget widget(file_name);
    widget.show();
#else
    QLabel note("OpenGL Support required");
    note.show();
#endif
    return app.exec();
  }
  catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Try `" << argv[0] << " --help' for more information."
              << std::endl;
    return -1;
  }
  return 0;
}
