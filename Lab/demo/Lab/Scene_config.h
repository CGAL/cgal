#ifndef SCENE_CONFIG_H
#define SCENE_CONFIG_H

#include <QtCore/qglobal.h>

#ifdef demo_framework_EXPORTS
#  define scene_EXPORTS
#endif

#ifdef scene_EXPORTS
#  define SCENE_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_EXPORT Q_DECL_IMPORT
#endif

#endif // SCENE_CONFIG_H
