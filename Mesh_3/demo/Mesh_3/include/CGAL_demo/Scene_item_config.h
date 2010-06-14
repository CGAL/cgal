#ifndef SCENE_ITEM_CONFIG_H
#define SCENE_ITEM_CONFIG_H

#include <QtCore/qglobal.h>

#ifdef scene_item_EXPORTS
#  define SCENE_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_ITEM_EXPORT Q_DECL_IMPORT
#endif

#endif // SCENE_ITEM_CONFIG_H
