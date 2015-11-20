#ifndef MAINWINDOW_CONFIG_H
#define MAINWINDOW_CONFIG_H

#include <QtCore/qglobal.h>

#ifdef polyhedron_demo_EXPORTS
#  define mainwindow_EXPORTS
#endif

#ifdef mainwindow_EXPORTS
#  define MAINWINDOW_EXPORT Q_DECL_EXPORT
#else
#  define MAINWINDOW_EXPORT Q_DECL_IMPORT
#endif

#endif // MAINWINDOW_CONFIG_H
