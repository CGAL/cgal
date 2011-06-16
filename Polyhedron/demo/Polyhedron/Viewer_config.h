#ifndef VIEWER_CONFIG_H
#define VIEWER_CONFIG_H

#include <QtCore/qglobal.h>

#ifdef demo_framework_EXPORTS
#  define viewer_EXPORTS
#endif

#ifdef viewer_EXPORTS
#  define VIEWER_EXPORT Q_DECL_EXPORT
#else
#  define VIEWER_EXPORT Q_DECL_IMPORT
#endif

#endif // VIEWER_CONFIG_H
