// Copyright (c) 2012-2015  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent RINEAU, Maxime Gimeno

#ifndef TEXTRENDERER_H
#define TEXTRENDERER_H

#include <CGAL/license/Three.h>


#include <QObject>
#include <QVector3D>

#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Scene_interface.h>

//!This class holds the properties of each line of text to be rendered.
class QVector3D;
namespace CGAL{
namespace Three{
class Scene_item;
}
}
class  VIEWER_EXPORT TextItem{
public :
  /*!
   * \brief The default constructor.
   * This is an overloaded functtion.
   */
    TextItem() {}
    /*!
     * \brief The construtor for the TextItem
     * \param p_x the X coordinate of the displayed text.
     * \param p_y the Y coordinate of the displayed text.
     * \param p_z the Z coordinate of the displayed text.
     * \param p_text the text to render.
     * \param p_3D
     * If true : the renderer will convert the coordinates into world coordinates.
     * If false : the renderer will display the text in screen coordinates, and ignore p_z
     * \param font the font used for the rendering.
     * \param p_color the color of the text.
     */
    TextItem(float p_x, float p_y, float p_z, QString p_text, bool p_3D = true, QFont font = QFont(), QColor p_color = Qt::black, bool always_visible = false)
        :x(p_x), y(p_y), z(p_z),_3D(p_3D), _is_always_visible(always_visible), m_text(p_text), m_font(font), m_color(p_color)
    {
       QFontMetrics fm(m_font);
       _width = float(fm.width(m_text));
       _height = float(fm.height());
    }
    QString text()const {return m_text;}
    //!Returns the position of the center of the text, in world coordinates.
    QVector3D position(){return QVector3D(x,y,z);}
    float width(){return _width;}
    float height(){return _height;}
    QFont font(){return m_font;}
    QColor color() {return m_color;}
    bool is_3D()const{return _3D;}
    bool is_always_visible(){return _is_always_visible;}
private:
    float x;
    float y;
    float z;
    float _width;
    float _height;
    bool _3D;
    bool _is_always_visible;
    QString m_text;
    QFont m_font;
    QColor m_color;

};//end class TextItem

class  VIEWER_EXPORT TextListItem{
public:
    TextListItem(CGAL::Three::Scene_item* pItem)
    {
      _list = QList<TextItem*>();
      _item = pItem;
    }

     CGAL::Three::Scene_item* item()const {return _item;}
    QList<TextItem*> textList()const {return _list;}
    void append(TextItem* ti) {_list.append(ti);}
    void clear(){_list.clear();}
    bool isEmpty()const {return _list.empty();}
    std::size_t size()const{return _list.size();}
private:
    CGAL::Three::Scene_item* _item;
    QList<TextItem*> _list;

};
//!This class draws all the textItems.
/*!
  * Projects each textItem from the world coordinates to the Screen coordinates
  * and draws it.
 */
class VIEWER_EXPORT TextRenderer : public QObject{
  Q_OBJECT
public:
    TextRenderer() : max_textItems(30000)
    {
    }
    //!Draws all the TextItems
    void draw(CGAL::Three::Viewer_interface* viewer);
    //!Adds a TextItem to the local list.
    void addText(TextItem*);
    //!Adds a TextListItem to the global list.
    void addTextList(TextListItem*);
    //!Creates a new TextItem and adds it to the local list.
    void addText(float p_x, float p_y, float p_z, QString p_text, bool p_3D = true,  QFont font = QFont(), QColor p_color = Qt::black);
    //!Removes a TextItem from the local list.
    void removeText(TextItem*);
    //!Removes a TextItemList from the global list.
    void removeTextList(TextListItem*);
    //!Returns the local list of TextItems. This is the renderer's default list.
    QList<TextItem*> getLocalTextItems(){return local_textItems;}
    //!Returns the global list of TextItems. This is the list that is fed by pre-filled lists of TextItems (such as global Polyhedron IDs).
    QList<TextListItem*> items() const{return textItems;}
    //!Gives the renderer a Scene, needed to determine which Ids must be drawn.
    void setScene(CGAL::Three::Scene_interface* p_scene){scene = p_scene;}
    int getMax_textItems()const{return max_textItems;}
    void setMax(int max){max_textItems = max;}

Q_SIGNALS:
    void sendMessage(QString message, int ms_delay =2000);
private:
    QList<TextListItem*> textItems;
    CGAL::Three::Scene_interface *scene;
    QList<TextItem*> local_textItems;
    int max_textItems;

};//end class TextRenderer
#endif // TEXTRENDERER_H
