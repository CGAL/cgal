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

class QVector3D;
namespace CGAL{
namespace Three{
class Scene_item;
}
}
//!A TextItem is a string with properties, like coordinates, font and color, that is
//! to be rendered in the viewer.
class  VIEWER_EXPORT TextItem{
public :
  /*!
   * \brief The default constructor.
   */
    TextItem() {}
    /*!
     * \brief The construtor for the TextItem
     * \param p_x, p_y, p_z the coordinates of the TextItem.
     * \param p_text the text to render.
     * \param p_3D
     * If true : the TextRenderer will convert the coordinates into world coordinates.
     * If false : the TextRenderer will display the text in screen coordinates, and ignore p_z
     * \param font the font used for the rendering.
     * \param p_color the color of the text.
     * \param always_visible overrides Viewer_interface::testDisplayId() if true;
     */
    TextItem(float p_x, float p_y, float p_z, QString p_text, bool p_3D = true, QFont font = QFont(), QColor p_color = Qt::black, bool always_visible = false)
        :x(p_x), y(p_y), z(p_z),_3D(p_3D), _is_always_visible(always_visible), m_text(p_text), m_font(font), m_color(p_color)
    {
       QFontMetrics fm(m_font);
       _width = float(fm.width(m_text));
       _height = float(fm.height());
    }
    //!\brief Accessor for the string
    //!
    //! @returns the item's string
    QString text()const {return m_text;}
    //!\brief The center of the item
    //!
    //!@returns the position of the center of the text, in world coordinates.
    QVector3D position(){return QVector3D(x,y,z);}
    //!The text's width in pixels
    float width(){return _width;}
    //!The text's height in pixels
    float height(){return _height;}
    //!Accessor for the text's font
    QFont font(){return m_font;}
    //!Accessor for the text's color
    QColor color() {return m_color;}
    //!Specifies if the item's coordinates are World or Screen coordinates
    //! @returns true if they are World coordinates
    //! @returns false if they are Screen coordinates
    bool is_3D()const{return _3D;}
    //!Specifies if the item may be hidden
    //! @returns true if it may not
    //! @returns false if it may
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

//! This is a list of TextItems associated with a CGAL::Three::Scene_item.
class  VIEWER_EXPORT TextListItem{
public:
  //!\brief The Constructor
  //! \param pItem the CGAL::Three::Scene_item associated to the list.
    TextListItem(CGAL::Three::Scene_item* pItem)
    {
      _list = QList<TextItem*>();
      _item = pItem;
    }

    //!The associated CGAL::Three::Scene_item.
    CGAL::Three::Scene_item* item()const {return _item;}
    //!The list of TextItems
    QList<TextItem*> textList()const {return _list;}
    //!Adds a `ti` to the list.
    void append(TextItem* ti) {_list.append(ti);}
    //!Clears the list of all the TextItems.
    void clear(){_list.clear();}
    //!Tells if the list there is any TextItem in the list.
    //! \returns true if there is no TextItem in the list
    //! \returns false if there is any TextItem in the list
    bool isEmpty()const {return _list.empty();}
    //!The size of the list
    std::size_t size()const{return _list.size();}
private:
    CGAL::Three::Scene_item* _item;
    QList<TextItem*> _list;

};
/*!
  * The TextRender projects each TextItem from the world coordinates to the Screen coordinates
  * and draws it.
 */
class VIEWER_EXPORT TextRenderer : public QObject{
  Q_OBJECT
public:
    TextRenderer() : max_textItems(10000)
    {
    }
    //!Draws all the `TextItem`s
    void draw(CGAL::Three::Viewer_interface* viewer);
    //!\brief Adds a single TextItem to TextRenderer::local_textItems
    //!
    //! @see addText(float p_x, float p_y, float p_z, QString p_text, bool p_3D = true,  QFont font = QFont(), QColor p_color = Qt::black)
    void addText(TextItem*);
    //!\brief Creates a new TextItem in TextRenderer::local_textItems
    //!
    //!This is a version of addText(TextItem*) that creates the TextItem on the fly.
    //! @see addText(TextItem*)
    void addText(float p_x, float p_y, float p_z, QString p_text, bool p_3D = true,  QFont font = QFont(), QColor p_color = Qt::black);
    //!\brief Adds a TextListItem to TextRenderer::textItems
    //!
    void addTextList(TextListItem*);
    //!Removes a `textItem` from TextRenderer::local_textItems
    //! @attention the memory of the TextItem is not de-allocated.
    //! @see addText(TextItem*)
    void removeText(TextItem* textItem);
    //!Removes a TextItemList from TextRenderer::textItems
    //! //! @attention the memory of the TextItems is not de-allocated.
    void removeTextList(TextListItem*);
    //!The local TextListItem.
    QList<TextItem*> getLocalTextItems(){return local_textItems;}
    //!The global list of `TextListItem`s.
    QList<TextListItem*> items() const{return textItems;}
    //!Gives the renderer a Scene.
    void setScene(CGAL::Three::Scene_interface* p_scene){scene = p_scene;}
    //! The maximum `TextItem`s that can be displayed at once.
    //! @see setMax()
    int getMax_textItems()const{return max_textItems;}
    //!Sets the maximum `TextItem`s that can be displayed at once.
    //! @see getMax_textItems()
    void setMax(int max){max_textItems = max;}

Q_SIGNALS:
    //!Emit this to print `message` on the viewer for `ms_delay` ms.
    void sendMessage(QString message, int ms_delay =2000);
protected:
    //!\brief list of `TextListItem`s
    //!
    //! It holds lists of correlated `TextItem`s that are associated to a
    //! CGAL::Three::Scene_item, like primitiveIds.
    QList<TextListItem*> textItems;
    //!\brief List of `TextItem`s
    //!
    //! Usually fed by the viewer, it holds the text informations from the
    //! viewer that are displayed directly on the screen, like the fps,
    //! the distances, etc.
    QList<TextItem*> local_textItems;
private:
    CGAL::Three::Scene_interface *scene;
    int max_textItems;

};//end class TextRenderer
#endif // TEXTRENDERER_H
