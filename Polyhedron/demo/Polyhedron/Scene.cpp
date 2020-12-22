#include "config.h"
#include "Scene.h"

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_print_item_interface.h>
#include <CGAL/Three/Scene_transparent_interface.h>
#include <CGAL/Three/Scene_zoomable_item_interface.h>
#include <CGAL/Three/Three.h>

#include  <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_print_item_interface.h>
#include <CGAL/Three/Viewer_interface.h>


#include <QObject>
#include <QMetaObject>
#include <QString>
#include <QGLWidget>
#include <QEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QColorDialog>
#include <QApplication>
#include <QPointer>
#include <QList>
#include <QAbstractProxyModel>
#include <QMimeData>
#include <QOpenGLFramebufferObject>


Scene::Scene(QObject* parent)
    : QStandardItemModel(parent),
      selected_item(-1),
      item_A(-1),
      item_B(-1)
{

    connect(this, SIGNAL(selectionRay(double, double, double,
                                      double, double, double)),
            this, SLOT(setSelectionRay(double, double, double,
                                       double, double, double)));
    connect(this, SIGNAL(indexErased(Scene_interface::Item_id)),
              this, SLOT(adjustIds(Scene_interface::Item_id)));
    picked = false;
    gl_init = false;
    dont_emit_changes = false;

}
Scene::Item_id
Scene::addItem(CGAL::Three::Scene_item* item)
{
    Bbox bbox_before = bbox();
    m_entries.push_back(item);
    Item_id id = m_entries.size() - 1;
    connect(item, SIGNAL(itemChanged()),
            this, SLOT(itemChanged()));
    connect(item, SIGNAL(itemVisibilityChanged()),
            this, SLOT(itemVisibilityChanged()));
    connect(item, SIGNAL(redraw()),
            this, SLOT(callDraw()));
    if(item->isFinite()
            && !item->isEmpty()
            && bbox_before + item->bbox() != bbox_before
            )
    {
        Q_EMIT updated_bbox(true);
    }
    QList<QStandardItem*> list;
    for(int i=0; i<5; i++)
    {
        list<<new QStandardItem();
        list.at(i)->setEditable(false);
    }
    invisibleRootItem()->appendRow(list);
    for(int i=0; i<5; i++){
        index_map[list.at(i)->index()] = m_entries.size() -1;
    }
    Q_EMIT updated();
    children.push_back(id);
    Q_EMIT newItem(id);
    CGAL::Three::Scene_group_item* group =
            qobject_cast<CGAL::Three::Scene_group_item*>(item);
    if(group)
        addGroup(group);
    //init the item for the mainViewer to avoid using unexisting
    //VAOs if the mainViewer is not the first to be drawn.
    QOpenGLFramebufferObject* fbo = CGAL::Three::Three::mainViewer()->depthPeelingFbo();
    CGAL::Three::Three::mainViewer()->setDepthPeelingFbo(nullptr);//to prevent crashing as the fbo is not initialized in this call.
    item->draw(CGAL::Three::Three::mainViewer());
    item->drawEdges(CGAL::Three::Three::mainViewer());
    item->drawPoints(CGAL::Three::Three::mainViewer());
    CGAL::Three::Three::mainViewer()->setDepthPeelingFbo(fbo);
    if(group)
       m_groups.append(id);
    return id;
}

CGAL::Three::Scene_item*
Scene::replaceItem(Scene::Item_id index, CGAL::Three::Scene_item* item, bool emit_item_about_to_be_destroyed)
{
    if(index < 0 || index >= m_entries.size())
        return 0;

    connect(item, SIGNAL(itemChanged()),
            this, SLOT(itemChanged()));
    connect(item, SIGNAL(itemVisibilityChanged()),
            this, SLOT(itemVisibilityChanged()));
    connect(item, SIGNAL(redraw()),
            this, SLOT(callDraw()));
    CGAL::Three::Scene_group_item* group =
            qobject_cast<CGAL::Three::Scene_group_item*>(m_entries[index]);
    QList<Scene_item*> group_children;
    if(group)
    {
      m_groups.removeAll(index);
      Q_FOREACH(Item_id id, group->getChildren())
      {
        CGAL::Three::Scene_item* child = group->getChild(id);
        group->unlockChild(child);
        group_children << child;
      }
    }
    CGAL::Three::Scene_group_item* parent = m_entries[index]->parentGroup();
    bool is_locked = false;
    if(parent)
    {
      is_locked = parent->isChildLocked(m_entries[index]);
      parent->unlockChild(m_entries[index]);
      parent->removeChild(m_entries[index]);
    }
    std::swap(m_entries[index], item);
    if(parent)
    {
      changeGroup(m_entries[index], parent);
      if(is_locked)
        parent->lockChild(m_entries[index]);
    }

    Q_EMIT newItem(index);
    if ( item->isFinite() && !item->isEmpty() &&
         m_entries[index]->isFinite() && !m_entries[index]->isEmpty() &&
         item->bbox()!=m_entries[index]->bbox() )
    {
      Q_EMIT updated_bbox(true);
    }

    if(emit_item_about_to_be_destroyed) {
      Q_EMIT itemAboutToBeDestroyed(item);
      item->aboutToBeDestroyed();
    }

    Q_EMIT updated();
    group =
            qobject_cast<CGAL::Three::Scene_group_item*>(m_entries[index]);
    if(group)
    {
        addGroup(group);
        m_groups.append(index);
    }
    itemChanged(index);
    Q_EMIT restoreCollapsedState();
    redraw_model();
    Q_EMIT selectionChanged(index);
    Q_FOREACH(Scene_item* child, group_children)
    {
      erase(item_id(child));
    }
    return item;
}

Scene::Item_id
Scene::erase(Scene::Item_id index)
{
  CGAL::Three::Scene_item* item = m_entries[index];
  if(qobject_cast<Scene_group_item*>(item))
  {
    setSelectedItemsList(QList<Scene_interface::Item_id>()<<item_id(item));
    return erase(selectionIndices());
  }
  m_groups.removeAll(index);
  if(item->parentGroup()
     && item->parentGroup()->isChildLocked(item))
    return -1;
  //clears the Scene_view
  clear();
  index_map.clear();
  if(index < 0 || index >= m_entries.size())
    return -1;
  if(item->parentGroup())
    item->parentGroup()->removeChild(item);

  //removes the item from all groups that contain it
  Item_id removed_item = item_id(item);
  children.removeAll(removed_item);
  indexErased(removed_item);
    m_entries.removeAll(item);
  Q_EMIT itemAboutToBeDestroyed(item);
  item->aboutToBeDestroyed();
  item->deleteLater();
  selected_item = -1;
  //re-creates the Scene_view
  Q_FOREACH(Item_id id, children)
  {
    organize_items(this->item(id), invisibleRootItem(), 0);
  }
  QStandardItemModel::beginResetModel();
  Q_EMIT updated();
  QStandardItemModel::endResetModel();
  Q_EMIT restoreCollapsedState();
  if(--index >= 0)
    return index;
  if(!m_entries.isEmpty())
    return 0;
  return -1;
}

int
Scene::erase(QList<int> indices)
{
  if(indices.empty())
    return -1;
  QList<CGAL::Three::Scene_item*> to_be_removed;
  int max_index = -1;
  Q_FOREACH(int index, indices) {
    if(index < 0 || index >= m_entries.size())
      continue;

    max_index = (std::max)(max_index, index);
    CGAL::Three::Scene_item* item = m_entries[index];
    if(item->parentGroup()
       && item->parentGroup()->isChildLocked(item))
      if(!indices.contains(item_id(item->parentGroup())))
        continue;
    Scene_group_item* group = qobject_cast<Scene_group_item*>(item);
    if(group)
    {
      Q_FOREACH(Item_id id, group->getChildren())
      {
        CGAL::Three::Scene_item* child = group->getChild(id);
        if(!to_be_removed.contains(child))
          to_be_removed.push_back(child);
      }
    }
    if(!to_be_removed.contains(item))
      to_be_removed.push_back(item);
  }

  Q_FOREACH(Scene_item* item, to_be_removed) {
    Item_id removed_item = item_id(item);
    if(removed_item == -1) //case of the selection_item, for example.
      continue;
    if(item->parentGroup())
      item->parentGroup()->removeChild(item);
    children.removeAll(removed_item);
    indexErased(removed_item);
    m_groups.removeAll(removed_item);
    m_entries.removeAll(item);

    Q_EMIT itemAboutToBeDestroyed(item);
    item->aboutToBeDestroyed();
    item->deleteLater();
  }
  clear();
  index_map.clear();
  selected_item = -1;
  Q_FOREACH(Item_id id, children)
  {
    organize_items(item(id), invisibleRootItem(), 0);
  }
  QStandardItemModel::beginResetModel();
  Q_EMIT updated();
  QStandardItemModel::endResetModel();
  Q_EMIT restoreCollapsedState();

  int index = max_index + 1 - indices.size();
  if(index >= m_entries.size()) {
    index = m_entries.size() - 1;
  }
  if(index >= 0)
    return index;
  if(!m_entries.isEmpty())
    return 0;
  return -1;

}

void Scene::remove_item_from_groups(Scene_item* item)
{
    CGAL::Three::Scene_group_item* group = item->parentGroup();
    if(group)
    {
        group->removeChild(item);
        children.push_back(item_id(item));
    }
}
Scene::~Scene()
{
  Q_FOREACH(CGAL::QGLViewer* viewer, CGAL::QGLViewer::QGLViewerPool())
  {
    removeViewer(static_cast<CGAL::Three::Viewer_interface*>(viewer));
    viewer->setProperty("is_destroyed", true);
  }
  Q_FOREACH(QOpenGLVertexArrayObject* vao, vaos.values())
  {
    vao->destroy();
    delete vao;
  }
  Q_FOREACH(CGAL::Three::Scene_item* item_ptr, m_entries)
  {
    item_ptr->deleteLater();
  }
  m_entries.clear();
}

CGAL::Three::Scene_item*
Scene::item(Item_id index) const
{
    return m_entries.value(index); // QList::value checks bounds
}

Scene::Item_id
Scene::item_id(CGAL::Three::Scene_item* scene_item) const
{
    return m_entries.indexOf(scene_item);
}

int
Scene::numberOfEntries() const
{
    return m_entries.size();
}

// Duplicate a scene item.
// Return the ID of the new item (-1 on error).
Scene::Item_id
Scene::duplicate(Item_id index)
{
    if(index < 0 || index >= m_entries.size())
        return -1;

    const CGAL::Three::Scene_item* item = m_entries[index];
    CGAL::Three::Scene_item* new_item = item->clone();
    if(new_item) {
        new_item->setName(tr("%1 (copy)").arg(item->name()));
        new_item->setColor(item->color());
        new_item->setVisible(item->visible());
        addItem(new_item);
        return m_entries.size() - 1;
    }
    else
        return -1;
}

void Scene::initializeGL(CGAL::Three::Viewer_interface* viewer)
{

  //Vertex source code
  const char vertex_source[] =
  {
    "#version 150                                 \n"
    "in vec4 vertex;                \n"
    "in vec2 v_texCoord;            \n"
    "uniform mat4 projection_matrix;       \n"
    "out vec2 f_texCoord;              \n"
    "void main(void)                             \n"
    "{                                           \n"
    "  f_texCoord = v_texCoord;                  \n"
    "  gl_Position = projection_matrix * vertex; \n"
    "}                                           \n"

  };

  const char vertex_source_comp[] =
  {
    "attribute highp vec4 vertex;                \n"
    "attribute highp vec2 v_texCoord;            \n"
    "uniform highp mat4 projection_matrix;       \n"
    "varying highp vec2 f_texCoord;              \n"
    "void main(void)                             \n"
    "{                                           \n"
    "  f_texCoord = v_texCoord;                  \n"
    "  gl_Position = projection_matrix * vertex; \n"
    "}                                           \n"

  };
  //Fragment source code
  const char fragment_source[] =
  {
    "#version 150                                                            \n"
    "in vec2 f_texCoord;                                         \n"
    "out vec4 out_color ; \n"
    "uniform sampler2D s_texture;                                             \n"
    "void main(void)                                                        \n"
    "{                                                                      \n"
    "  out_color = texture(s_texture, f_texCoord); \n"
    "}                                                                      \n"
  };
  const char fragment_source_comp[] =
  {
    "varying highp vec2 f_texCoord;                                         \n"
    "uniform sampler2D texture;                                             \n"
    "void main(void)                                                        \n"
    "{                                                                      \n"
    "  gl_FragColor = texture2D(texture, f_texCoord); \n"
    "}                                                                      \n"
  };


  QOpenGLShader vertex_shader(QOpenGLShader::Vertex);
  QOpenGLShader fragment_shader(QOpenGLShader::Fragment);
  if(viewer->isOpenGL_4_3())
  {
    if(!vertex_shader.compileSourceCode(vertex_source))
    {
      std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }

    if(!fragment_shader.compileSourceCode(fragment_source))
    {
      std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }
  }
  else
  {
    if(!vertex_shader.compileSourceCode(vertex_source_comp))
    {
      std::cerr<<"Compiling vertex source FAILED"<<std::endl;
    }

    if(!fragment_shader.compileSourceCode(fragment_source_comp))
    {
      std::cerr<<"Compiling fragmentsource FAILED"<<std::endl;
    }
  }

  if(!program.addShader(&vertex_shader))
  {
    std::cerr<<"adding vertex shader FAILED"<<std::endl;
  }
  if(!program.addShader(&fragment_shader))
  {
    std::cerr<<"adding fragment shader FAILED"<<std::endl;
  }
  if(!program.link())
  {
    //std::cerr<<"linking Program FAILED"<<std::endl;
    qDebug() << program.log();
  }
  points[0] = -1.0f; points[1] = -1.0f; points[2] = 0.0f;
  points[3] = 1.0f; points[4] = 1.0f; points[5] = 0.0f;
  points[6] = 1.0f; points[7] = -1.0f; points[8] = 0.0f;
  points[9] = -1.0f; points[10] = -1.0f; points[11] = 0.0f;
  points[12] = -1.0f; points[13] = 1.0f; points[14] = 0.0f;
  points[15] = 1.0f; points[16] = 1.0f; points[17] = 0.0f;

  uvs[0] = 0.0f; uvs[1] = 0.0f;
  uvs[2] = 1.0f; uvs[3] = 1.0f;
  uvs[4] = 1.0f; uvs[5] = 0.0f;
  uvs[6] = 0.0f; uvs[7] = 0.0f;
  uvs[8] = 0.0f; uvs[9] = 1.0f;
  uvs[10] = 1.0f; uvs[11] = 1.0f;

  vbo[0].create();
  vbo[1].create();

  viewer->makeCurrent();
  vaos[viewer] = new QOpenGLVertexArrayObject();
  vaos[viewer]->create();
  program.bind();
  vaos[viewer]->bind();
  vbo[0].bind();
  vbo[0].allocate(points, 18 * sizeof(float));
  program.enableAttributeArray("vertex");
  program.setAttributeArray("vertex", GL_FLOAT, 0, 3);
  vbo[0].release();

  vbo[1].bind();
  vbo[1].allocate(uvs, 12 * sizeof(float));
  program.enableAttributeArray("v_texCoord");
  program.setAttributeArray("v_texCoord", GL_FLOAT, 0, 2);
  vbo[1].release();
  vaos[viewer]->release();
  program.release();
  gl_init = true;
}

void Scene::s_itemAboutToBeDestroyed(CGAL::Three::Scene_item *rmv_itm)
{
 Q_FOREACH(CGAL::Three::Scene_item* item, m_entries)
 {
   if(item == rmv_itm)
     item->itemAboutToBeDestroyed(item);
 }
}
bool
Scene::keyPressEvent(QKeyEvent* e){
    bool res=false;
    for (QList<int>::iterator it=selected_items_list.begin(),endit=selected_items_list.end();
         it!=endit;++it)
    {
        CGAL::Three::Scene_item* item=m_entries[*it];
        res |= item->keyPressEvent(e);
    }
    return res;
}

void
Scene::draw(CGAL::Three::Viewer_interface* viewer)
{
    draw_aux(false, viewer);
}
void
Scene::drawWithNames(CGAL::Three::Viewer_interface* viewer)
{
    draw_aux(true, viewer);
}

bool item_should_be_skipped_in_draw(Scene_item* item) {
  if(!item->visible()) return true;
  if(item->has_group == 0) return false;
  Scene_group_item* group = item->parentGroup();
  while(group != 0) {
    if(!group->visible()) return false;
    group = group->parentGroup();
  }
  return true;
}


void Scene::renderScene(const QList<Scene_interface::Item_id> &items,
                        Viewer_interface *viewer,
                        QMap<float, int>& picked_item_IDs,
                        bool with_names,
                        int pass,
                        bool writing_depth,
                        QOpenGLFramebufferObject *fbo)
{
  viewer->setCurrentPass(pass);
  viewer->setDepthWriting(writing_depth);
  viewer->setDepthPeelingFbo(fbo);
  Q_FOREACH(Scene_interface::Item_id index, items)
  {
    CGAL::Three::Scene_item& item = *m_entries[index];
    CGAL::Three::Scene_group_item* group =
        qobject_cast<CGAL::Three::Scene_group_item*>(&item);
    if(index == selected_item || selected_items_list.contains(index))
    {
      item.selection_changed(true);
    }
    else
    {
      item.selection_changed(false);
    }
    if(group ||item.visible())
    {
      if( group || item.renderingMode() == Flat || item.renderingMode() == FlatPlusEdges || item.renderingMode() == Gouraud || item.renderingMode() == GouraudPlusEdges )
      {
        if(with_names) {
          viewer->glClearDepthf(1.0);
          viewer->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        }
        item.draw(viewer);
        if(with_names) {

          //    read depth buffer at pick location;
          float depth = 1.0;
          viewer->glReadPixels(picked_pixel.x(),viewer->camera()->screenHeight()-1-picked_pixel.y(),1,1,GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
          if (depth != 1.0)
          {
            //add object to list of picked objects;
            picked_item_IDs[depth] = index;
          }
        }
      }
      if(group)
        group->renderChildren(viewer, picked_item_IDs, picked_pixel, with_names);
    }
  }
}

void Scene::renderWireScene(const QList<Scene_interface::Item_id> &items,
                            Viewer_interface *viewer,
                            QMap<float, int>& picked_item_IDs,
                            bool with_names)
{
  Q_FOREACH(Scene_interface::Item_id index, items)
   {
     CGAL::Three::Scene_item& item = *m_entries[index];
     CGAL::Three::Scene_group_item* group =
         qobject_cast<CGAL::Three::Scene_group_item*>(&item);
     if(index == selected_item || selected_items_list.contains(index))
     {
         item.selection_changed(true);
     }
     else
     {
         item.selection_changed(false);
     }

     if(group ||item.visible())
     {
       if( group || (!with_names && item.renderingMode() == FlatPlusEdges )
          || item.renderingMode() == Wireframe
          || item.renderingMode() == PointsPlusNormals
          || item.renderingMode() == GouraudPlusEdges)
       {
         if(with_names) {
           viewer->glClearDepthf(1.0);
           viewer->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
         }
         viewer->setGlPointSize(2.f);
         item.drawEdges(viewer);
       }
       else{
           if( item.renderingMode() == PointsPlusNormals ){
               viewer->setGlPointSize(2.f);
               if(index == selected_item || selected_items_list.contains(index))
               {

                 item.selection_changed(true);
               }
               else
               {

                 item.selection_changed(false);
               }
               item.drawEdges(viewer);
           }
       }

       if((item.renderingMode() == Wireframe || item.renderingMode() == PointsPlusNormals )
          && with_names)
       {

         //    read depth buffer at pick location;
         float depth = 1.0;
         viewer->glReadPixels(picked_pixel.x(),viewer->camera()->screenHeight()-1-picked_pixel.y(),1,1,GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
         if (depth != 1.0)
         {
           //add object to list of picked objects;
           picked_item_IDs[depth] = index;
         }
       }
     }
   }
}

void Scene::renderPointScene(const QList<Scene_interface::Item_id> &items,
                             Viewer_interface *viewer,
                             QMap<float, int>& picked_item_IDs,
                             bool with_names)
{
  Q_FOREACH(Scene_interface::Item_id index, items)
  {
    CGAL::Three::Scene_item& item = *m_entries[index];
    CGAL::Three::Scene_group_item* group =
        qobject_cast<CGAL::Three::Scene_group_item*>(&item);
    if(group ||item.visible())
    {
      if(item.renderingMode() == Points && with_names) {
          viewer->glClearDepthf(1.0);
          viewer->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      }

      if(group || item.renderingMode() == Points  ||
         (item.renderingMode() == PointsPlusNormals)  ||
         (item.renderingMode() == ShadedPoints))
      {
        if(with_names) {
          viewer->glClearDepthf(1.0);
          viewer->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        }
        viewer->setGlPointSize(3.0f);
        item.drawPoints(viewer);
      }
      if(item.renderingMode() == Points && with_names) {
        //    read depth buffer at pick location;
        float depth = 1.0;
        viewer->glReadPixels(picked_pixel.x(),viewer->camera()->screenHeight()-1-picked_pixel.y(),1,1,GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
        if (depth != 1.0)
        {
          //add object to list of picked objects;
          picked_item_IDs[depth] = index;
        }
      }
    }
  }
}


 bool Scene::has_alpha()
 {
   Q_FOREACH(Scene_item* item, m_entries)
     if(item->alpha() != 1.0f)
       return true;
   return false;
 }
void
Scene::draw_aux(bool with_names, CGAL::Three::Viewer_interface* viewer)
{
    QMap<float, int> picked_item_IDs;
    if(with_names)
      viewer->glEnable(GL_DEPTH_TEST);
    if(!gl_init)
        initializeGL(viewer);
    //treat opaque items first to ensure that when two items are the same, but only one is opaque,
    //the item stays opaque
    QList<Item_id> opaque_items;
    QList<Item_id> transparent_items;
    Q_FOREACH(Item_id id, children)
    {
      Scene_item* item = m_entries[id];
      Scene_group_item* group = qobject_cast<Scene_group_item*>(item);
      bool is_transparent=false;
      if(item->alpha() != 1.0f)
        is_transparent = true;
      else if(group)
      {
        for(const auto& child : group->getChildren())
        {
          if(group->getChild(child)->alpha() < 1.0f)
          {
            is_transparent = true;
            break;
          }
        }
      }
      if(!is_transparent)
        opaque_items.push_back(id);
      else
        transparent_items.push_back(id);
    }
    renderScene(children, viewer, picked_item_IDs, with_names, -1, false, NULL);
    if(with_names)
    {
      //here we get the selected point, before erasing the depth buffer. We store it
      //in a dynamic property as a QList<double>. If there is some alpha, the
      //depth buffer is altered, and the picking will return true even when it is
      // performed in the background, when it should return false. To avoid that,
      // we distinguish the case were there is no alpha, to let the viewer
      //perform it, and the case where the pixel is not found. In the first case,
      //we erase the property, in the latter we return an empty list.
      //According ot that, in the viewer, either we perform the picking, either we do nothing.
      if(has_alpha()) {
        bool found = false;
        CGAL::qglviewer::Vec point = viewer->camera()->pointUnderPixel(picked_pixel, found) - viewer->offset();
        if(found){
          QList<QVariant> picked_point;
          picked_point <<point.x
                      <<point.y
                     <<point.z;
          viewer->setProperty("picked_point", picked_point);
        }
        else{
          viewer->setProperty("picked_point", QList<QVariant>());
        }
      }
      else {
        viewer->setProperty("picked_point", {});
      }
    }
    if(!with_names && has_alpha())
    {
      std::vector<QOpenGLFramebufferObject*> fbos;
      std::vector<QOpenGLFramebufferObject*> depth_test;
      QColor background = viewer->backgroundColor();
      fbos.resize((int)viewer->total_pass());
      depth_test.resize((int)viewer->total_pass()-1);

      //first pass
      fbos[0] = new QOpenGLFramebufferObject(viewer->width(), viewer->height(),QOpenGLFramebufferObject::Depth, GL_TEXTURE_2D, GL_RGBA32F);
      fbos[0]->bind();
      viewer->glDisable(GL_BLEND);
      viewer->glEnable(GL_DEPTH_TEST);
      viewer->glDepthFunc(GL_LESS);
      viewer->glClearColor(0.0f,
                           0.0f,
                           0.0f,
                           0.0f);
      viewer->glClearDepthf(1);
      viewer->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      renderScene(opaque_items, viewer, picked_item_IDs, false, 0,false, NULL);
      renderScene(transparent_items, viewer, picked_item_IDs, false, 0,false, NULL);
      fbos[0]->release();
      depth_test[0] = new QOpenGLFramebufferObject(viewer->width(), viewer->height(),QOpenGLFramebufferObject::Depth, GL_TEXTURE_2D, GL_RGBA32F);
      depth_test[0]->bind();
      viewer->glDisable(GL_BLEND);
      viewer->glEnable(GL_DEPTH_TEST);
      viewer->glDepthFunc(GL_LESS);
      viewer->glClearColor(0.0f,
                           0.0f,
                           0.0f,
                           0.0f);
      viewer->glClearDepthf(1);
      viewer->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      renderScene(opaque_items, viewer, picked_item_IDs, false, 0,true, NULL);
      renderScene(transparent_items, viewer, picked_item_IDs, false, 0,true, NULL);
      depth_test[0]->release();

      //other passes
      for(int i=1; i<viewer->total_pass()-1; ++i)
      {
        fbos[i] = new QOpenGLFramebufferObject(viewer->width(), viewer->height(),QOpenGLFramebufferObject::Depth, GL_TEXTURE_2D, GL_RGBA32F);
        fbos[i]->bind();
        viewer->glDisable(GL_BLEND);
        viewer->glEnable(GL_DEPTH_TEST);
        viewer->glDepthFunc(GL_LESS);
        viewer->glClearColor(0.0f,
                             0.0f,
                             0.0f,
                             0.0f);
        viewer->glClearDepthf(1);
        viewer->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        renderWireScene(children, viewer, picked_item_IDs, false);
        renderPointScene(children, viewer, picked_item_IDs, false);
        renderScene(opaque_items     , viewer, picked_item_IDs, false, i, false, depth_test[i-1]);
        renderScene(transparent_items, viewer, picked_item_IDs, false, i, false, depth_test[i-1]);
        fbos[i]->release();

        depth_test[i] = new QOpenGLFramebufferObject(viewer->width(), viewer->height(),QOpenGLFramebufferObject::Depth, GL_TEXTURE_2D, GL_RGBA32F);
        depth_test[i]->bind();
        viewer->glDisable(GL_BLEND);
        viewer->glEnable(GL_DEPTH_TEST);
        viewer->glDepthFunc(GL_LESS);
        viewer->glClearColor(0.0f,
                             0.0f,
                             0.0f,
                             0.0f);
        viewer->glClearDepthf(1);
        viewer->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        renderScene(opaque_items     , viewer, picked_item_IDs, false, i, true, depth_test[i-1]);
        renderScene(transparent_items, viewer, picked_item_IDs, false, i, true, depth_test[i-1]);
        depth_test[i]->release();
      }


      //last pass
      fbos[(int)viewer->total_pass()-1] = new QOpenGLFramebufferObject(viewer->width(), viewer->height(),QOpenGLFramebufferObject::Depth, GL_TEXTURE_2D, GL_RGBA32F);
      fbos[(int)viewer->total_pass()-1]->bind();
      viewer->glDisable(GL_BLEND);
      viewer->glEnable(GL_DEPTH_TEST);
      viewer->glDepthFunc(GL_LESS);
      viewer->glClearColor(0.0f,
                           0.0f,
                           0.0f,
                           0.0f);
      viewer->glClearDepthf(1);
      viewer->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      renderScene(opaque_items     , viewer, picked_item_IDs, false, (int)viewer->total_pass()-1, false, depth_test[(int)viewer->total_pass()-2]);
      renderScene(transparent_items, viewer, picked_item_IDs, false, (int)viewer->total_pass()-1, false, depth_test[(int)viewer->total_pass()-2]);
      fbos[(int)viewer->total_pass()-1]->release();
      if(viewer->getStoredFrameBuffer() != NULL)
        viewer->getStoredFrameBuffer()->bind();

      //blending
      program.bind();
      vaos[viewer]->bind();
      viewer->glClearColor((GLclampf)background.redF(),
                           (GLclampf)background.greenF(),
                           (GLclampf)background.blueF(),
                           0.0f);
      viewer->glDisable(GL_DEPTH_TEST);
      viewer->glClear(GL_COLOR_BUFFER_BIT);
      viewer->glEnable(GL_BLEND);
      viewer->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

      QMatrix4x4 proj_mat;
      proj_mat.setToIdentity();
      proj_mat.ortho(-1,1,-1,1,0,1);
      program.setUniformValue("projection_matrix", proj_mat);
      for(int i=0; i< (int)viewer->total_pass()-1; ++i)
        delete depth_test[i];
      for(int i = (int)viewer->total_pass()-1; i>=0; --i)
      {
        viewer->glBindTexture(GL_TEXTURE_2D, fbos[i]->texture());
        viewer->glDrawArrays(GL_TRIANGLES,0,static_cast<GLsizei>(6));
        delete fbos[i];
      }
      viewer->glDisable(GL_BLEND);
      viewer->glEnable(GL_DEPTH_TEST);
      vaos[viewer]->release();
      program.release();
    }

    viewer->glDepthFunc(GL_LEQUAL);
    // Wireframe OpenGL drawing
    renderWireScene(children, viewer, picked_item_IDs, with_names);
    // Points OpenGL drawing
    renderPointScene(children, viewer, picked_item_IDs, with_names);

    if(with_names)
    {
        QList<float> depths = picked_item_IDs.keys();
        if(!depths.isEmpty())
        {
            std::sort(depths.begin(), depths.end());
            int id = picked_item_IDs[depths.first()];
            setSelectedItemIndex(id);
            viewer->setSelectedName(id);

        }
    }
    if(with_names)
        picked = true;
    else
        picked = false;
    //scrolls the sceneView to the selected item's line.
    if(picked)
    {
        Q_EMIT(itemPicked(index_map.key(mainSelectionIndex())));
    }
    Q_EMIT drawFinished();
}

// workaround for Qt-4.2 (see above)
#undef lighter
QVariant
Scene::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
    {
        return QVariant();
    }

    int id = index_map[index];
    if(id < 0 || id >= m_entries.size())
        return QVariant();
    if(role == ::Qt::ToolTipRole)
    {
        return m_entries[id]->toolTip();
    }
    switch(index.column())
    {
    case ColorColumn:
        if(role == ::Qt::DecorationRole)
            return m_entries.value(id)->color();
        break;
    case NameColumn:
        if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
            return m_entries.value(id)->name();
        if(role == ::Qt::FontRole)
            return m_entries.value(id)->font();
        break;
    case RenderingModeColumn:
        if(role == ::Qt::DisplayRole) {
            return m_entries.value(id)->renderingModeName();
        }
        else if(role == ::Qt::EditRole) {
            return static_cast<int>(m_entries.value(id)->renderingMode());
        }
        else if(role == ::Qt::TextAlignmentRole) {
            return ::Qt::AlignCenter;
        }
        break;
    case ABColumn:
        if(role == ::Qt::DisplayRole) {
            if(id == item_A)
                return "A";
            if(id == item_B)
                return "B";
        }
        else if(role == ::Qt::TextAlignmentRole) {
            return ::Qt::AlignLeft;
        }
        break;
    case VisibleColumn:
        if(role == ::Qt::DisplayRole || role == ::Qt::EditRole)
            return m_entries.value(id)->visible();
        break;
    default:
        return QVariant();
    }
    return QVariant();
}

QVariant
Scene::headerData ( int section, ::Qt::Orientation orientation, int role ) const
{
    if(orientation == ::Qt::Horizontal)  {
        if (role == ::Qt::DisplayRole)
        {
            switch(section)
            {
            case NameColumn:
                return tr("Name");
                break;
            case ColorColumn:
                return tr("#");
                break;
            case RenderingModeColumn:
                return tr("Mode");
            case ABColumn:
                return tr("A/B");
                break;
            case VisibleColumn:
                return tr("View");
                break;
            default:
                return QVariant();
            }
        }
        else if(role == ::Qt::ToolTipRole) {
            if(section == RenderingModeColumn) {
                return tr("Rendering mode (points/wireframe/flat/flat+edges/Gouraud)");
            }
            else if(section == ABColumn) {
                return tr("Selection A/Selection B");
            }
        }
    }
    return QStandardItemModel::headerData(section, orientation, role);
}

Qt::ItemFlags
Scene::flags ( const QModelIndex & index ) const
{
    if (index.isValid() && index.column() == NameColumn) {
        return QStandardItemModel::flags(index) | ::Qt::ItemIsEditable;
    }
    else {
        return QStandardItemModel::flags(index);
    }
}

bool
Scene::setData(const QModelIndex &index,
               const QVariant &value,
               int role)
{

    if( role != ::Qt::EditRole || !index.isValid() )
        return false;

    int id = index_map[index];
    if(id < 0 || id >= m_entries.size()){
        return false;
    }

    CGAL::Three::Scene_item* item = m_entries[id];

    if(!item) return false;
    switch(index.column())
    {
    case NameColumn:
        item->setName(value.toString());
    Q_EMIT dataChanged(index, index);
        return true;
        break;
    case ColorColumn:
      if(selectionIndices().contains(item_id(item)))
        Q_FOREACH(Item_id item_index, selectionIndices())
          this->item(item_index)->setColor(value.value<QColor>());
      else
        item->setColor(value.value<QColor>());
    Q_EMIT dataChanged(index, index);
        return true;
        break;
    case RenderingModeColumn:
    {
        RenderingMode rendering_mode = static_cast<RenderingMode>(value.toInt());
        // Find next supported rendering mode
        int counter = 0;
        while ( ! item->supportsRenderingMode(rendering_mode)
                )
        {
            rendering_mode = static_cast<RenderingMode>( (rendering_mode+1) % NumberOfRenderingMode );
            if(counter++ == NumberOfRenderingMode)
              break;
        }
        item->setRenderingMode(rendering_mode);
        QModelIndex nindex = createIndex(m_entries.size()-1,RenderingModeColumn+1);
    Q_EMIT dataChanged(index, nindex);
        return true;
        break;
    }
    case VisibleColumn:
        item->setVisible(value.toBool());
    Q_EMIT dataChanged(index, createIndex(m_entries.size()-1,VisibleColumn+1));
        return true;
    default:
        return false;
    }
    return false;
}

bool Scene::dropMimeData(const QMimeData * /*data*/,
                         Qt::DropAction /*action*/,
                         int /*row*/,
                         int /*column*/,
                         const QModelIndex &parent)
{
    //gets the moving items
    QList<Scene_item*> items;
    QList<int> groups_children;

    //get IDs of all children of selected groups
    Q_FOREACH(int i, selected_items_list)
    {
      CGAL::Three::Scene_group_item* group =
          qobject_cast<CGAL::Three::Scene_group_item*>(item(i));
      if(group)
        Q_FOREACH(Item_id id, group->getChildren())
        {
          CGAL::Three::Scene_item* child = item(id);
          groups_children << item_id(child);
        }
    }
    // Insure that children of selected groups will not be added twice
    Q_FOREACH(int i, selected_items_list)
    {
      if(!groups_children.contains(i))
      {
        items << item(i);
      }
    }
    //Gets the group at the drop position
    CGAL::Three::Scene_group_item* group = NULL;
    if(parent.isValid())
        group = qobject_cast<CGAL::Three::Scene_group_item*>(this->item(index_map[parent]));
    bool one_contained = false;
    if(group)
    {
      Q_FOREACH(int id, selected_items_list)
        if(group->getChildren().contains(id))
        {
          one_contained = true;
          break;

        }
    }
    //if the drop item is not a group_item or if it already contains the item, then the drop action must be ignored
    if(!group ||one_contained)
    {
      //unless the drop zone is empty, which means the item should be removed from all groups.
      if(!parent.isValid())
      {
        Q_FOREACH(Scene_item* item, items)
        {
          if(item->parentGroup())
          {
            item->parentGroup()->removeChild(item);
            addChild(item);
          }
        }
        redraw_model();
        return true;
      }
      return false;
    }
    Q_FOREACH(Scene_item* item, items)
      changeGroup(item, group);
    redraw_model();
    return true;
}

//todo : if a group is selected, don't treat it's children.
bool Scene::sort_lists(QVector<QList<int> >&sorted_lists, bool up)
{
  QVector<int> group_found;
  Q_FOREACH(int i, selectionIndices())
  {
    Scene_item* item = this->item(i);
    if(item->has_group == 0)
    {
      sorted_lists.first().push_back(i);
    }
    else
    {
      int group_id = item_id(item->parentGroup());
      if(group_found.contains(group_id))
        sorted_lists[group_id].push_back(i);
      else
      {
        group_found.push_back(group_id);
        if(sorted_lists.size() < group_id+1)
          sorted_lists.resize(group_id+1);
        sorted_lists[group_id].push_back(i);
      }
    }
  }
  //iterate the first list to find the groups that are selected and remove the corresponding
  //sub lists.
  //todo: do that for each group. (treat subgroups)
  for(int i = 0; i< sorted_lists.first().size(); ++i)
  {
    Scene_group_item* group = qobject_cast<Scene_group_item*>(this->item(sorted_lists.first()[i]));
    if(group && ! group->getChildren().isEmpty())
    {
      sorted_lists[sorted_lists.first()[i]].clear();
    }
  }
  std::sort(sorted_lists.first().begin(), sorted_lists.first().end(),
            [this](int a, int b) {
    return children.indexOf(a) < children.indexOf(b);
});
  if(!sorted_lists.first().isEmpty())
  {
    if(up &&  children.indexOf(sorted_lists.first().first()) == 0)
      return false;
    else if(!up &&  children.indexOf(sorted_lists.first().last()) == children.size() -1)
      return false;
  }
  for(int i=1; i<sorted_lists.size(); ++i)
  {
    QList<int>& list = sorted_lists[i];
    if(list.isEmpty())
      continue;
    Scene_group_item* group = qobject_cast<Scene_group_item*>(this->item(i));
    if(!group)
      continue;
    std::sort(list.begin(), list.end(),
              [group](int a, int b) {
      return group->getChildren().indexOf(a) < group->getChildren().indexOf(b);
  });
    if(up && group->getChildren().indexOf(list.first()) == 0)
      return false;
    else if(!up && group->getChildren().indexOf(list.last()) == group->getChildren().size()-1)
      return false;
  }
  return true;
}
void Scene::moveRowUp()
{
  if(selectionIndices().isEmpty())
    return;
  QVector<QList<int> >sorted_lists(1);
  QList<int> to_select;
  //sort lists according to the indices of each item in its container (scene or group)
  //if moving one up would put it out of range, then we stop and do nothing.
  if(!sort_lists(sorted_lists, true))
    return;

  for(int i=0; i<sorted_lists.first().size(); ++i)
  {
    Item_id selected_id = sorted_lists.first()[i];
    Scene_item* selected_item = item(selected_id);
    if(!selected_item)
      return;
    if(index_map.key(selected_id).row() > 0)
    {
      //if not in group
      QModelIndex baseId = index_map.key(selected_id);
      int newId = children.indexOf(
            index_map.value(index(baseId.row()-1, baseId.column(),baseId.parent()))) ;
      children.move(children.indexOf(selected_id), newId);
      redraw_model();
      to_select.append(m_entries.indexOf(selected_item));
    }
  }
  for(int i=1; i<sorted_lists.size(); ++i)
  {
    for(int j = 0; j< sorted_lists[i].size(); ++j)
    {
      Item_id selected_id = sorted_lists[i][j];
      Scene_item* selected_item = item(selected_id);
      if(!selected_item)
        return;
      if(index_map.key(selected_id).row() > 0)
      {
        Scene_group_item* group = selected_item->parentGroup();
        if(group)
        {
          int id = group->getChildren().indexOf(item_id(selected_item));
          group->moveUp(id);
          redraw_model();
          to_select.append(m_entries.indexOf(selected_item));
        }
      }
    }
  }
  if(!to_select.isEmpty()){
    selectionChanged(to_select);
  }
}
void Scene::moveRowDown()
{
  if(selectionIndices().isEmpty())
    return;
  QVector<QList<int> >sorted_lists(1);
  QList<int> to_select;
  //sort lists according to the indices of each item in its container (scene or group)
  //if moving one up would put it out of range, then we stop and do nothing.
  if(!sort_lists(sorted_lists, false))
    return;
  for(int i=sorted_lists.first().size()-1; i>=0; --i)
  {
    Item_id selected_id = sorted_lists.first()[i];
    Scene_item* selected_item = item(selected_id);
    if(!selected_item)
      return;
    if(index_map.key(selected_id).row() < rowCount(index_map.key(selected_id).parent())-1)
    {
        //if not in group
        QModelIndex baseId = index_map.key(selected_id);
        int newId = children.indexOf(
              index_map.value(index(baseId.row()+1, baseId.column(),baseId.parent()))) ;
        children.move(children.indexOf(selected_id), newId);

      redraw_model();
      to_select.prepend(m_entries.indexOf(selected_item));
    }
  }
  for(int i=1; i<sorted_lists.size(); ++i){
    if(sorted_lists[i].isEmpty())
      continue;
    for(int j = sorted_lists[i].size()-1; j >=0; --j)
    {
      Item_id selected_id = sorted_lists[i][j];
      Scene_item* selected_item = item(selected_id);
      if(!selected_item)
        return;
      if(index_map.key(selected_id).row() < rowCount(index_map.key(selected_id).parent())-1)
      {
        if(item(selected_id)->has_group >0)
        {
          Scene_group_item* group = selected_item->parentGroup();
          if(group)
          {
            int id = group->getChildren().indexOf(item_id(selected_item));
            group->moveDown(id);
          }
        }
        redraw_model();
        to_select.prepend(m_entries.indexOf(selected_item));
      }
    }
  }
  if(!to_select.isEmpty()){
    selectionChanged(to_select);
  }
}
Scene::Item_id Scene::mainSelectionIndex() const {
    return (selectionIndices().size() == 1) ? selected_item : -1;
}

QList<int> Scene::selectionIndices() const {
    return selected_items_list;
}

int Scene::selectionAindex() const {
    return item_A;
}

int Scene::selectionBindex() const {
    return item_B;
}

QItemSelection Scene::createSelection(int i)
{
    return QItemSelection(index_map.keys(i).at(0),
                          index_map.keys(i).at(4));
}

QItemSelection Scene::createSelection(QList<int> is)
{
    QItemSelection sel;
    Q_FOREACH(int i, is)
      sel.select(index_map.keys(i).at(0),
                 index_map.keys(i).at(4));
    return sel;
}

QItemSelection Scene::createSelectionAll()
{
  //it is not possible to directly create a selection with items that have different parents, so
  //we do it iteratively.
  QItemSelection sel;
  sel.select(this->createIndex(0, 0),
             this->createIndex(m_entries.size(), LastColumn));
  for(const auto& gid : m_groups)
  {
    CGAL::Three::Scene_group_item* group =
        qobject_cast<CGAL::Three::Scene_group_item*>(item(gid));
    sel.select(index_map.keys(group->getChildren().first()).at(0),
               index_map.keys(group->getChildren().last()).at(4));
  }
  return sel;
}

void Scene::itemChanged()
{
    CGAL::Three::Scene_item* item = qobject_cast<CGAL::Three::Scene_item*>(sender());
    if(item)
        itemChanged(item);
}

void Scene::itemChanged(Item_id i)
{
  if(dont_emit_changes)
    return;
  if(i < 0 || i >= m_entries.size())
    return;

  Q_EMIT dataChanged(this->createIndex(i, 0),
                     this->createIndex(i, LastColumn));
}

void Scene::itemChanged(CGAL::Three::Scene_item*item )
{
  if(dont_emit_changes)
    return;
  itemChanged(item_id(item));
}

void Scene::allItemsChanged()
{
  Q_EMIT dataChanged(this->createIndex(0, 0),
                     this->createIndex(m_entries.size() - 1, LastColumn));
}

void Scene::itemVisibilityChanged()
{
    CGAL::Three::Scene_item* item = qobject_cast<CGAL::Three::Scene_item*>(sender());
    if(item)
        itemVisibilityChanged(item);
}

void Scene::itemVisibilityChanged(CGAL::Three::Scene_item* item)
{
  if(item->isFinite()
     && !item->isEmpty())
  {
    //does not recenter
    if(visibility_recentering_enabled){
      Q_EMIT updated_bbox(true);

    }
  }
}


bool SceneDelegate::editorEvent(QEvent *event, QAbstractItemModel *model,
                                const QStyleOptionViewItem &option,
                                const QModelIndex &index)
{
    QAbstractProxyModel* proxyModel = dynamic_cast<QAbstractProxyModel*>(model);
    Q_ASSERT(proxyModel);
    Scene *scene = dynamic_cast<Scene*>(proxyModel->sourceModel());
    Q_ASSERT(scene);
    int id = scene->index_map[proxyModel->mapToSource(index)];
    switch(index.column()) {
    case Scene::VisibleColumn:
        if (event->type() == QEvent::MouseButtonPress) {
            QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
            if(mouseEvent->button() == ::Qt::LeftButton) {
                int x = mouseEvent->pos().x() - option.rect.x();
                if(x >= (option.rect.width() - size)/2 &&
                        x <= (option.rect.width() + size)/2) {
                    model->setData(index, !model->data(index).toBool());
                }
            }
            return false; //so that the selection can change
        }
        return true;
        break;
    case Scene::ColorColumn:
        if (event->type() == QEvent::MouseButtonPress) {
            QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
            if(mouseEvent->button() == ::Qt::LeftButton) {
                QColor color =
                        QColorDialog::getColor(model->data(index).value<QColor>(),
                                               0/*,
                                               tr("Select color"),
                                               QColorDialog::ShowAlphaChannel*/);
                if (color.isValid()) {
                    model->setData(index, color );
                }
            }
        }
        else if(event->type() == QEvent::MouseButtonDblClick) {
            return true; // block double-click
        }
        return false;
        break;
    case Scene::RenderingModeColumn:
        if (event->type() == QEvent::MouseButtonPress) {
            QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
            if(mouseEvent->button() == ::Qt::LeftButton) {
                // Switch rendering mode
                /*RenderingMode*/int rendering_mode = model->data(index, ::Qt::EditRole).toInt();
                rendering_mode = (rendering_mode+1) % NumberOfRenderingMode;
                model->setData(index, rendering_mode);
            }
        }
        else if(event->type() == QEvent::MouseButtonDblClick) {
            return true; // block double-click
        }
        return false;
        break;
    case Scene::ABColumn:
        if (event->type() == QEvent::MouseButtonPress) {
            if(id == scene->item_B) {
                scene->item_A = id;
                scene->item_B = -1;
            }
            else if(id == scene->item_A) {
                scene->item_B = id;
                scene->item_A = -1;
            }
            else if(scene->item_A == -1) {
                scene->item_A = id;
            }
            else {
                scene->item_B = id;
            }
            scene->dataChanged(scene->createIndex(0, Scene::ABColumn),
                               scene->createIndex(scene->rowCount() - 1, Scene::ABColumn));
        }
        return false;
        break;
    default:
        return QItemDelegate::editorEvent(event, model, option, index);
    }
}

void SceneDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option,
                          const QModelIndex &index) const
{
    QModelIndex test = proxy->mapToSource(index);
    if (index.column() != Scene::VisibleColumn) {
        QItemDelegate::paint(painter, option, index);
    } else {
        const QAbstractItemModel *model = index.model();

        QPalette::ColorGroup cg = (option.state & QStyle::State_Enabled) ?
                    (option.state & QStyle::State_Active) ? QPalette::Normal : QPalette::Inactive : QPalette::Disabled;

        if (option.state & QStyle::State_Selected)
            painter->fillRect(option.rect, option.palette.color(cg, QPalette::Highlight));
        bool checked = model->data(index, ::Qt::DisplayRole).toBool();
        int width = option.rect.width();
        int height = option.rect.height();
        size = (std::min)(width, height);
        int x = option.rect.x() + (option.rect.width() / 2) - (size / 2);;
        int y = option.rect.y() + (option.rect.height() / 2) - (size / 2);
        if(test.row()>=0 && test.row()<scene->m_entries.size()){

            if(checked) {
                painter->drawPixmap(x, y, checkOnPixmap.scaled(QSize(size, size),
                                                               ::Qt::KeepAspectRatio,
                                                               ::Qt::SmoothTransformation));
            }
            else {
                painter->drawPixmap(x, y, checkOffPixmap.scaled(QSize(size, size),
                                                                ::Qt::KeepAspectRatio,
                                                                ::Qt::SmoothTransformation));
            }
        }
        drawFocus(painter, option, option.rect); // since we draw the grid ourselves
    }
}

void Scene::setItemVisible(int index, bool b)
{
    if( index < 0 || index >= m_entries.size() )
        return;
    m_entries[index]->setVisible(b);
  Q_EMIT dataChanged(this->createIndex(index, VisibleColumn),
                     this->createIndex(index, VisibleColumn));
}

void Scene::setSelectionRay(double orig_x,
                            double orig_y,
                            double orig_z,
                            double dir_x,
                            double dir_y,
                            double dir_z)
{
    CGAL::Three::Scene_item* item = this->item(selected_item);
    if(item) item->select(orig_x,
                          orig_y,
                          orig_z,
                          dir_x,
                          dir_y,
                          dir_z);
}

void Scene::setItemA(int i)
{
    item_A = i;
    if(item_A == item_B)
    {
        item_B = -1;
    }
  Q_EMIT dataChanged(this->createIndex(0, ABColumn),
                     this->createIndex(m_entries.size()-1, ABColumn));
}

void Scene::setItemB(int i)
{
    item_B = i;
    if(item_A == item_B)
    {
        item_A = -1;
    }
  Q_EMIT updated();
  Q_EMIT dataChanged(this->createIndex(0, ABColumn),
                     this->createIndex(m_entries.size()-1, ABColumn));
}

Scene::Bbox Scene::bbox() const
{
    if(m_entries.empty())
        return Bbox(0,0,0,0,0,0);

    bool bbox_initialized = false;
    Bbox bbox = Bbox(0,0,0,0,0,0);
    Q_FOREACH(CGAL::Three::Scene_item* item, m_entries)
    {
        if(item->isFinite() && !item->isEmpty() && item->visible()) {
            if(bbox_initialized) {

                bbox = bbox + item->bbox();
            }
            else {
                bbox = item->bbox();
                bbox_initialized = true;

            }
        }

    }
    return bbox;
}

QList<Scene_item*> Scene::item_entries() const
{
    return m_entries;
}
void Scene::redraw_model()
{
    //makes the hierarchy in the tree
    //clears the model
    clear();
    index_map.clear();
    //fills the model
    Q_FOREACH(Item_id id, children)
    {
        organize_items(m_entries[id], invisibleRootItem(), 0);
    }
    Q_EMIT restoreCollapsedState();
}
void Scene::changeGroup(Scene_item *item, CGAL::Three::Scene_group_item *target_group)
{
    //remove item from the containing group if any
    if(item->parentGroup())
    {
      if(item->parentGroup()->isChildLocked(item))
        return;
      item->parentGroup()->removeChild(item);
      children.push_back(item_id(item));
    }
      else
      {
        children.removeAll(item_id(item));
      }
    //add the item to the target group
    target_group->addChild(item);
    item->moveToGroup(target_group);
    redraw_model();
    Q_EMIT updated();
}

void Scene::printPrimitiveId(QPoint point, CGAL::Three::Viewer_interface* viewer)
{
  Scene_item *it = item(mainSelectionIndex());
  if(it)
  {
    //Only call printPrimitiveId if the item is a Scene_print_item_interface
    Scene_print_item_interface* item= qobject_cast<Scene_print_item_interface*>(it);
    if(item)
      item->printPrimitiveId(point, viewer);
  }
}
void Scene::printVertexIds()
{
  Scene_item *it = item(mainSelectionIndex());
  if(it)
  {
    Scene_print_item_interface* item= qobject_cast<Scene_print_item_interface*>(it);
    if(item)
      item->printVertexIds();
  }
}

void Scene::printEdgeIds()
{
  Scene_item *it = item(mainSelectionIndex());
  if(it)
  {
    //Only call printEdgeIds if the item is a Scene_print_item_interface
    Scene_print_item_interface* item= qobject_cast<Scene_print_item_interface*>(it);
    if(item)
      item->printEdgeIds();
  }
}

void Scene::printFaceIds()
{
  Scene_item *it = item(mainSelectionIndex());
  if(it)
  {
    //Only call printFaceIds if the item is a Scene_print_item_interface
    Scene_print_item_interface* item= qobject_cast<Scene_print_item_interface*>(it);
    if(item)
      item->printFaceIds();
  }
}

void Scene::printAllIds()
{
  Scene_item *it = item(mainSelectionIndex());
  if(it)
  {
    //Only call printFaceIds if the item is a Scene_print_item_interface
    Scene_print_item_interface* item= qobject_cast<Scene_print_item_interface*>(it);
    if(item)
      item->printAllIds();
  }
}
void Scene::updatePrimitiveIds(CGAL::Three::Scene_item* it)
{
  if(it)
  {
    Scene_print_item_interface* item= qobject_cast<Scene_print_item_interface*>(it);
    if(item)
    {
      //As this function works as a toggle, the first call hides the ids and the second one shows them again,
      //thereby triggering their re-computation.
      item->printVertexIds();
      item->printVertexIds();

      item->printEdgeIds();
      item->printEdgeIds();

      item->printFaceIds();
      item->printFaceIds();
    }
  }
}
bool Scene::testDisplayId(double x, double y, double z, CGAL::Three::Viewer_interface* viewer)
{
    CGAL::Three::Scene_item *i = item(mainSelectionIndex());
    if(!i)
      return false;
    Scene_print_item_interface* spit= qobject_cast<Scene_print_item_interface*>(i);
    if(spit && i->visible())
    {
        bool res = spit->testDisplayId(x,y,z, viewer);
        return res;
    }
    else
      return false;
}
#include "Scene_find_items.h"

void Scene::organize_items(Scene_item* item, QStandardItem* root, int loop)
{
    if(item->has_group <= loop)
    {
        QList<QStandardItem*> list;
        for(int i=0; i<5; i++)
        {
            list<<new QStandardItem();
            list.at(i)->setEditable(false);

        }
        root->appendRow(list);
        for(int i=0; i<5; i++){
            index_map[list.at(i)->index()] = m_entries.indexOf(item);
        }
        CGAL::Three::Scene_group_item* group =
                qobject_cast<CGAL::Three::Scene_group_item*>(item);
        if(group)
        {
          Q_FOREACH(Item_id id, group->getChildren())
          {
            CGAL::Three::Scene_item* child = group->getChild(id);
                organize_items(child, list.first(), loop+1);
            }
        }
    }
}

void Scene::setExpanded(QModelIndex id)
{
    CGAL::Three::Scene_group_item* group =
            qobject_cast<CGAL::Three::Scene_group_item*>(item(getIdFromModelIndex(id)));
    if(group)
    {
        group->setExpanded(true);
    }
}
void Scene::setCollapsed(QModelIndex id)
{
    CGAL::Three::Scene_group_item* group =
            qobject_cast<CGAL::Three::Scene_group_item*>(item(getIdFromModelIndex(id)));
    if(group)
    {
        group->setExpanded(false);
    }
}

int Scene::getIdFromModelIndex(QModelIndex modelId)const
{
    return index_map.value(modelId);
}

QList<QModelIndex> Scene::getModelIndexFromId(int id) const
{
    return index_map.keys(id);
}

void Scene::addGroup(Scene_group_item* group)
{
    connect(this, SIGNAL(drawFinished()), group, SLOT(resetDraw()));
    connect(this, SIGNAL(indexErased(Scene_interface::Item_id)),
                group, SLOT(adjustIds(Scene_interface::Item_id)));
}

namespace scene { namespace details {

Q_DECL_EXPORT
CGAL::Three::Scene_item*
findItem(const CGAL::Three::Scene_interface* scene_interface,
         const QMetaObject& metaobj,
         QString name, Scene_item_name_fn_ptr fn) {
    const Scene* scene = dynamic_cast<const Scene*>(scene_interface);
    if(!scene) return 0;
    Q_FOREACH(CGAL::Three::Scene_item* item, scene->entries()) {
       CGAL::Three::Scene_item* ptr = qobject_cast<CGAL::Three::Scene_item*>(metaobj.cast(item));
        if(ptr && ((ptr->*fn)() == name)) return ptr;
    }
    return 0;
}

Q_DECL_EXPORT
QList<CGAL::Three::Scene_item*>
findItems(const CGAL::Three::Scene_interface* scene_interface,

          const QMetaObject&,
          QString name, Scene_item_name_fn_ptr fn)
{
    const Scene* scene = dynamic_cast<const Scene*>(scene_interface);
    QList<CGAL::Three::Scene_item*> list;
    if(!scene) return list;

    Q_FOREACH(CGAL::Three::Scene_item* item, scene->entries()) {
        CGAL::Three::Scene_item* ptr = qobject_cast<CGAL::Three::Scene_item*>(item);
        if(ptr && ((ptr->*fn)() == name)) {
            list << ptr;
        }
    }
    return list;
}

} // end namespace details
                } // end namespace scene

void Scene::zoomToPosition(QPoint point, Viewer_interface *viewer)
{
  for(int i=0; i<numberOfEntries(); ++i)
  {
    if(!item(i)->visible())
      continue;
    Scene_zoomable_item_interface* zoom_item = qobject_cast<Scene_zoomable_item_interface*>(item(i));
    if(zoom_item)
    {
      zoom_item->zoomToPosition(point, viewer);
    }
  }
}

void Scene::adjustIds(Item_id removed_id)
{
  for(int i = 0; i < children.size(); ++i)
  {
    if(children[i] >= removed_id)
      --children[i];
  }
  for(int i = removed_id; i < numberOfEntries(); ++i)
  {
    m_entries[i]->setId(i-1);//the signal is emitted before m_entries is amputed from the item, so new id is current id -1.
  }
}

void Scene::computeBbox()
{
  if(m_entries.empty())
  {
    last_bbox = Bbox(0,0,0,0,0,0);
    return;
  }

  bool bbox_initialized = false;
  Bbox bbox = Bbox(0,0,0,0,0,0);
  Q_FOREACH(CGAL::Three::Scene_item* item, m_entries)
  {
    if(item->isFinite() && !item->isEmpty() ) {
      if(bbox_initialized) {

        bbox = bbox + item->bbox();
      }
      else {
        bbox = item->bbox();
        bbox_initialized = true;

      }
    }

  }
  last_bbox = bbox;
}

void Scene::newViewer(Viewer_interface *viewer)
{
  initGL(viewer);
  Q_FOREACH(Scene_item* item, m_entries)
  {
    item->newViewer(viewer);
  }
}

void Scene::removeViewer(Viewer_interface *viewer)
{
 //already destroyed;
  if(viewer->property("is_destroyed").toBool())
    return;

  vaos[viewer]->destroy();
  vaos[viewer]->deleteLater();
  vaos.remove(viewer);
  Q_FOREACH(Scene_item* item, m_entries)
  {
    item->removeViewer(viewer);
  }
}

void Scene::initGL(Viewer_interface *viewer)
{
  viewer->makeCurrent();
  vaos[viewer] = new QOpenGLVertexArrayObject();
  vaos[viewer]->create();
  program.bind();
  vaos[viewer]->bind();
  vbo[0].bind();
  vbo[0].allocate(points, 18 * sizeof(float));
  program.enableAttributeArray("vertex");
  program.setAttributeArray("vertex", GL_FLOAT, 0, 3);
  vbo[0].release();

  vbo[1].bind();
  vbo[1].allocate(uvs, 12 * sizeof(float));
  program.enableAttributeArray("v_texCoord");
  program.setAttributeArray("v_texCoord", GL_FLOAT, 0, 2);
  vbo[1].release();
  vaos[viewer]->release();
  program.release();
}

void Scene::callDraw(){
  Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
  {
    qobject_cast<Viewer_interface*>(v)->update();
  }
}

void Scene::enableVisibilityRecentering(bool b)
{
  visibility_recentering_enabled = b;
}

void Scene::addChild(Scene_item *item)
{
  children.push_back(item_id(item));
}
