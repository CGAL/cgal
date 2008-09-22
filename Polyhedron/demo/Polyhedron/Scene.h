#ifndef SCENE_H
#define SCENE_H

#include <QtOpenGL/qgl.h>
#include <QAbstractListModel>
#include <QString>
#include <QColor>
#include <QList>
#include <QItemDelegate>
#include <QPixmap>
#include <QItemSelection>

#include "Nef_type_fwd.h" // declares Nef_polyhedron
#include "Polyhedron_type_fwd.h" // declares Polyhedron
#include "Textured_polyhedron_type_fwd.h" // declares textured polyhedron

#include <iostream>
#include <cmath>
#include <boost/variant.hpp>

#include "texture.h"

class QEvent;
class QMouseEvent;

class Scene  :
  public QAbstractListModel
{
  Q_OBJECT

  friend class SceneDelegate;
public:
  enum RenderingMode { Fill = 0, 
                       Wireframe, 
                       LastRenderingMode = Wireframe,
                       NumberOfRenderingMode = Wireframe+1};

  enum Columns { NameColumn = 0, 
                 ColorColumn, 
                 RenderingModeColumn, 
                 ActivatedColumn,
                 ABColumn,
                 LastColumn = ABColumn,
                 NumberOfColumns = LastColumn + 1};


  Scene(QObject*  parent);
  ~Scene();

  void addPolyhedron(Polyhedron* p,
                     QString name,
                     QColor color = defaultColor,
                     bool activated = true,
                     RenderingMode mode = Fill);

  void addTexPolyhedron(Textured_polyhedron* p,
                     QString name,
                     QColor color = defaultColor,
                     bool activated = true,
                     RenderingMode mode = Fill);

  void addNefPolyhedron(Nef_polyhedron* p,
			QString name,
			QColor color = defaultColor,
			bool activated = true,
			RenderingMode mode = Fill);

  int open(QString);  // Returns the index of the new polyhedra (-1 if
                      // error)
  bool save(int,QString); // Returns true upon successful save

  int erase(int);     // Returns the index of the polyhedra just before the
                      // one that is erased, or just after. Returns -1 if
                      // the list is empty.

  int duplicate(int); // Returns the index of the new polyhedra

  // Accessors (getters)
  int numberOfPolyhedra() const;
  Polyhedron* polyhedron(int i) const;
  Nef_polyhedron* nefPolyhedron(int i) const;
  Textured_polyhedron* texPolyhedron(int i) const;
  
  enum Entry_type { POLYHEDRON_ENTRY = 0,
                    NEF_ENTRY = 1,
		    TEX_POLYHEDRON_ENTRY = 2};
  Entry_type polyhedronType(int) const;

  QColor polyhedronColor(int) const;
  QString polyhedronName(int) const;
  bool isPolyhedronActivated(int) const;
  RenderingMode polyhedronRenderingMode(int) const;
  int selectionAindex() const;
  int selectionBindex() const;

  // for backward compatibility
  Polyhedron* getPolyhedron(int i) { return polyhedron(i); }

  // initializeGL() is called by Viewer::initializeGL()
  void initializeGL();
  // draw() is called by Viewer::draw()
  void draw(bool with_names = false);

  struct Bbox {
    double xmin, ymin, zmin;
    double xmax, ymax, zmax;
    Bbox(const double _xmin,const double _ymin,const double _zmin,
         const double _xmax,const double _ymax,const double _zmax)
	 : xmin(_xmin), ymin(_ymin), zmin(_zmin),
	   xmax(_xmax), ymax(_ymax), zmax(_zmax)
    {
    }
    Bbox()
	 : xmin(0.0), ymin(0.0), zmin(0.0),
	   xmax(1.0), ymax(1.0), zmax(1.0)
    {
    }
  };

  double len_diagonal()
  {
    Bbox box = bbox();
    double dx = box.xmax - box.xmin;
    double dy = box.ymax - box.ymin;
    double dz = box.zmax - box.zmin;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  // defined in Scene_nef_and_polyhedron_operations.cpp
  Bbox bbox();

  // QAbstractItemModel functions
  int rowCount ( const QModelIndex & parent = QModelIndex() ) const;
  int columnCount ( const QModelIndex & parent = QModelIndex() ) const;
  QVariant data ( const QModelIndex & index, int role = ::Qt::DisplayRole ) const;
  QVariant headerData ( int section, ::Qt::Orientation orientation, int role = ::Qt::DisplayRole ) const;
  ::Qt::ItemFlags flags ( const QModelIndex & index ) const;
  bool setData(const QModelIndex &index, const QVariant &value, int role);

  // auxiliary public function for QMainWindow
  QItemSelection createSelection(int i);
public slots:
  void polyhedronChanged(int i);
  void polyhedronChanged(Polyhedron*);
  void setSelectedItem(int i )
  {
    selected_item = i;
  };

  void setViewEdges(bool b) 
  {
    viewEdges = b;
    emit updated();
  }
  // Accessors (setters)
  void setPolyhedronActivated(int, bool b);
  void setPolyhedronA(int i);
  void setPolyhedronB(int i);

signals:
  void updated_bbox();
  void updated();

private:
  // functions that need to know the type Polyhedron
  // defined in Scene_polyhedron_operations.cpp
  QString polyhedronToolTip(int index) const;
  Polyhedron* new_polyhedron();
  Polyhedron* copy_polyhedron(Polyhedron* poly);
  void destroy_polyhedron(Polyhedron*);
  bool load_polyhedron(Polyhedron* poly, std::istream& in); // return true
                                                            // iif the
                                                            // loading is ok.
  bool save_polyhedron(Polyhedron* poly, std::ostream& out); // return true
                                                             // iif the
                                                             // save is ok.

  // functions that need to know the type Nef_polyhedron
  // defined in Scene_nef_polyhedron_operations.cpp
  QString nefPolyhedronToolTip(int index) const;
  Nef_polyhedron* new_nef_polyhedron();
  Nef_polyhedron* copy_nef_polyhedron(Nef_polyhedron* poly);
  void destroy_nef_polyhedron(Nef_polyhedron*);

  // functions defined in Scene_Textured_polyhedron_operations.cpp
  QString texPolyhedronToolTip(int index) const;
  Textured_polyhedron* new_tex_polyhedron(); 
  Textured_polyhedron* copy_tex_polyhedron(Textured_polyhedron* poly);
  void destroy_tex_polyhedron(Textured_polyhedron* poly);


private:
  static const QColor defaultColor; // defined in Scene.cpp

  typedef boost::variant<Polyhedron*, Nef_polyhedron*, Textured_polyhedron*> Polyhedron_ptr;

  struct Polyhedron_entry {
    Polyhedron_entry()
            : rendering_mode(Fill),
              display_list_built(false) {};

    Polyhedron_ptr polyhedron_ptr;
    QString name;
    QColor color;
    bool activated;
    RenderingMode rendering_mode;

    // display list
    unsigned int display_list;
    unsigned int display_list_for_edges;
    bool display_list_built;
  };

  Polyhedron_ptr copy_polyhedron_ptr(Polyhedron_ptr);

  void addEntry(Polyhedron_ptr p,
		QString name,
		QColor color = defaultColor,
		bool activated = true,
		RenderingMode mode = Fill);
  void destroyEntry(Polyhedron_entry&);
  void destroy_entry_ptr(Polyhedron_ptr);

  void draw(Polyhedron_entry& entry); // draw one entry
  void gl_render_facets(Polyhedron_ptr);


  typedef QList<Polyhedron_entry> Polyhedra;
  Polyhedra polyhedra;
  int selected_item;
  int item_A;
  int item_B;
  bool viewEdges;
  Texture texture;
}; // end class Scene

class SceneDelegate : public QItemDelegate
{
public:
  SceneDelegate(QObject * parent = 0)
    : QItemDelegate(parent),
      checkOnPixmap(":/cgal/icons/check-on.png"),
      checkOffPixmap(":/cgal/icons/check-off.png")
  {
  }

  bool editorEvent(QEvent *event, QAbstractItemModel *model,
                   const QStyleOptionViewItem &option,
                   const QModelIndex &index);
  void paint(QPainter *painter, const QStyleOptionViewItem &option,
             const QModelIndex &index) const;

private:
  QPixmap checkOnPixmap;
  QPixmap checkOffPixmap;
  mutable int size;
}; // end class SceneDelegate

#endif // SCENE_H
