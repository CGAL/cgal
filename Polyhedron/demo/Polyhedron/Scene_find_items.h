#ifndef SCENE_FIND_ITEMS
#define SCENE_FIND_ITEMS

#include <QObject>
#include <QMetaObject>
#include "Scene_item.h" // required, to have &Scene_item::name
#include "Scene_config.h"

class Scene_interface;

namespace scene {
namespace details {

typedef QString (Scene_item ::*Scene_item_name_fn_ptr)() const;

// Declaration only (defined in Scene.cpp)
SCENE_EXPORT 
Scene_item* 
findItem(const Scene_interface* scene_interface,
         const QMetaObject& metaobj,
         QString name, Scene_item_name_fn_ptr fn); 

// Declaration only (defined in Scene.cpp)
SCENE_EXPORT
QList<Scene_item*> 
findItems(const Scene_interface* scene_interface, 
          const QMetaObject& metaobj,
          QString name, Scene_item_name_fn_ptr fn); // fwd declaration

template <typename T>
T findItem(const Scene_interface* scene, QString name,
           Scene_item_name_fn_ptr fn) 
{
  return 
    static_cast<T>(findItem(scene,
                            reinterpret_cast<T>(0)->staticMetaObject,
                            name, fn));
}

template <typename T>
QList<T> findItems(const Scene_interface* scene, QString name, 
                   Scene_item_name_fn_ptr fn)
{
  QList<Scene_item*> void_list = 
    findItems(scene, reinterpret_cast<T>(0)->staticMetaObject,
              name, fn);
  QList<T> list;
  Q_FOREACH(Scene_item* ptr, void_list) {
    list << qobject_cast<T>(ptr);
  }
  return list;
}

} // end namespace details

// Searches

/** Search the first item that can be cast to T (T must be a pointer
    type), and called "name". If "name" is omitted, all names are
    accepted.
*/      
template <typename T>
T findItem(const Scene_interface* scene, 
           QString item_name = QString())
{ return details::findItem<T>(scene, item_name, &Scene_item::name); }

/** Returns all items that can be cast to T (T must be a pointer
    type), and called "name". If "name" is omitted, all names are
    accepted.
*/      
template <typename T>
QList<T> findItems(const Scene_interface* scene, 
                   QString item_name = QString()) 
{ return details::findItems<T>(scene, item_name, &Scene_item::name); }

/** Search the first item that can be cast to T (T must be a pointer
    type), and that has objectName() equal to "name". If "name" is
    omitted, all names are accepted.
*/      
template <typename T>
T findItemByObjectName(const Scene_interface* scene, 
                       QString obj_name = QString()) 
{ return details::findItem<T>(scene, obj_name, &QObject::objectName); }

/** Returns all items that can be cast to T (T must be a pointer type),
    and have objectName() equal to "name". If "name" is omitted, all
    names are accepted.
*/      
template <typename T>
QList<T> findItemsByObjectName(const Scene_interface* scene, 
                               QString obj_name = QString()) 
{ return details::findItems<T>(scene, obj_name, &QObject::objectName); }


// template <typename T>
// T scene_findItem(const Scene* scene, QString name, Scene_item_name_fn_ptr fn) {
//   return 
//     static_cast<T>(scene_findItem(scene, name, fn, 
//                                   reinterpret_cast<T>(0)->staticMetaObject()));
// }

// template <typename T>
// QList<T> 
// scene_findItems(const Scene* scene, QString name, Scene_item_name_fn_ptr fn) {
//   QList<void*> void_list = 
//     scene_findItems(scene, name, fn, 
//                     reinterpret_cast<T>(0)->staticMetaObject());
//   QList<T> list;
//   Q_FOREACH(void* ptr, void_list) {
//     list << static_cast<T>(ptr);
//   }
//   return list;
// }

// } // end scene namespace 

} // end namespace scene



#endif
