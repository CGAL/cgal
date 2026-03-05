// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// Copyright (c) 2026  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s)     : Laurent Rineau


#include <CGAL/Bbox_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_voronoi_diagram_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_local_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/IO/File_poly.h>
#include <CGAL/IO/OBJ.h>
#include <CGAL/IO/write_VTU.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Mesh_2/Face_badness.h>
#include <CGAL/Object.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/DelaunayMeshTriangulationGraphicsItem.h>
#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/resources.h>
#include <CGAL/Qt/utility.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>

#include "DelaunayMeshInsertSeeds.h"
#include "TriangulationCircumcircle.h"

#include <QColor>
#include <QObject>
#include <QAction>
#include <QActionGroup>
#include <QApplication>
#include <QCheckBox>
#include <QColorDialog>
#include <QComboBox>
#include <QDialog>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFormLayout>
#include <QGraphicsScene>
#include <QGraphicsSceneHoverEvent>
#include <QGraphicsView>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QLabel>
#include <QLineEdit>
#include <QMainWindow>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QObject>
#include <QPainter>
#include <QPalette>
#include <QPushButton>
#include <QRegularExpression>
#include <QSettings>
#include <QSpinBox>
#include <QStatusBar>
#include <QString>
#include <QStringList>
#include <Qt>
#include <QTimer>
#include <QToolBar>
#include <QVBoxLayout>

#include <qgraphicsitem.h>
#include <qobjectdefs.h>
#include <qtmetamacros.h>

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <random>
#include <utility>
#include <vector>

using K1 = CGAL::Simple_cartesian<double>;
using Kernel = CGAL::Filtered_kernel<K1>;

struct K : public Kernel
{};

using FT = K::FT;
using Point_2 = K::Point_2;
using Segment_2 = K::Segment_2;
using Triangle_2 = K::Triangle_2;

using Vb = CGAL::Triangulation_vertex_base_2<K>;
using Fb = CGAL::Delaunay_mesh_face_base_2<K>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, Tds, CGAL::Exact_predicates_tag>;

using Criteria = CGAL::Delaunay_mesh_local_size_criteria_2<CDT>;
using Mesher = CGAL::Delaunay_mesher_2<CDT, Criteria>;

using Vertex_handle = CDT::Vertex_handle;
using Face_handle = CDT::Face_handle;
using Edge = CDT::Edge;

template <class CDTType> void read_constraints(CDTType& t, std::istream& f) {
  using Point = typename CDTType::Point;

  t.clear();

  int nedges = 0;
  f >> nedges;

  for(int n = 0; n < nedges; ++n) {
    Point p1, p2;
    f >> p1 >> p2;
    t.insert_constraint(p1, p2);
  }
}

template <class CDTType> void write_constraints(const CDTType& t, std::ostream& f) {
  int number_of_constrained_edges = 0;
  for(const auto& [face, edge_idx] : t.finite_edges()) {
    if(face->is_constrained(edge_idx)) {
      ++number_of_constrained_edges;
    }
  }

  f << number_of_constrained_edges << '\n';
  for(const auto& [face, edge_idx] : t.finite_edges()) {
    if(face->is_constrained(edge_idx)) {
      f << face->vertex(t.cw(edge_idx))->point() << " " << face->vertex(t.ccw(edge_idx))->point() << '\n';
    }
  }
}

// Helper range adaptors for Mesher API that only provides begin/end pairs
namespace {
auto bad_faces_range(const Mesher* m) {
  return CGAL::make_range(m->bad_faces_begin(), m->bad_faces_end());
}

auto encroached_edges_range(const Mesher* m) {
  return CGAL::make_range(m->encroached_edges_begin(), m->encroached_edges_end());
}

template <typename Clusters> auto clusters_vertices_range(const Clusters& c) {
  return CGAL::make_range(c.clusters_vertices_begin(), c.clusters_vertices_end());
}

QPen pen(QColor color,
         int width = 0,
         Qt::PenStyle style = Qt::SolidLine,
         Qt::PenCapStyle cap = Qt::RoundCap,
         Qt::PenJoinStyle join = Qt::RoundJoin) {
  QPen pen_{color, static_cast<qreal>(width), style, cap, join};
  pen_.setCosmetic(true);
  return pen_;
}

} // end anonymous namespace

class BadFacesGraphicsItem : public CGAL::Qt::GraphicsItem
{
public:
  BadFacesGraphicsItem(const CDT* cdt, const Mesher* mesher)
      : cdt_(cdt)
      , mesher_(mesher)
      , bounds_(-100.0, -100.0, 200.0, 200.0) {}

  void setMesher(const Mesher* mesher) {
    mesher_ = mesher;
    update();
  }

  void updateBounds(const QRectF& rect) {
    bounds_ = rect;
    prepareGeometryChange();
    update();
  }

  QRectF boundingRect() const override { return bounds_; }

  void paint(QPainter* painter, const QStyleOptionGraphicsItem*, QWidget*) override {
    if(mesher_ == nullptr || cdt_ == nullptr) {
      return;
    }

    CGAL::Qt::Converter<K> convert;

    painter->setBrush(brush_);
    painter->setPen(Qt::NoPen);

    for(auto fit : bad_faces_range(mesher_)) {
      Triangle_2 tri(fit->vertex(0)->point(), fit->vertex(1)->point(), fit->vertex(2)->point());
      painter->drawPolygon(convert(tri));
    }
  }

  void modelChanged() override {
    // do nothing
  }

  void setBrush(const QBrush& brush) {
    brush_ = brush;
    update();
  }

  QBrush brush() const { return brush_; }

private:
  const CDT* cdt_;
  const Mesher* mesher_;
  QRectF bounds_;
  QBrush brush_{QColor(0, 160, 0, 128)};
};

class EncroachedEdgesGraphicsItem : public CGAL::Qt::GraphicsItem
{
public:
  EncroachedEdgesGraphicsItem(const CDT* cdt, const Mesher* mesher)
      : cdt_(cdt)
      , mesher_(mesher)
      , bounds_(-100.0, -100.0, 200.0, 200.0) {}

  void setMesher(const Mesher* mesher) {
    mesher_ = mesher;
    update();
  }

  void updateBounds(const QRectF& rect) {
    bounds_ = rect;
    prepareGeometryChange();
    update();
  }

  QRectF boundingRect() const override { return bounds_; }

  void paint(QPainter* painter, const QStyleOptionGraphicsItem*, QWidget*) override {
    if(mesher_ == nullptr || cdt_ == nullptr) {
      return;
    }

    CGAL::Qt::Converter<K> convert;
    painter->setPen(pen_);

    for(const auto& e : encroached_edges_range(mesher_)) {
      Segment_2 s = cdt_->segment(e);
      painter->drawLine(convert(s));
    }
  }

  void modelChanged() override {
    // do nothing
  }

  void setPen(const QPen& pen) {
    pen_ = pen;
    update();
  }

  QPen pen() const { return pen_; }

private:
  const CDT* cdt_;
  const Mesher* mesher_;
  QRectF bounds_;
  QPen pen_{QColor(200, 0, 0), 0};
};

class ClustersGraphicsItem : public CGAL::Qt::GraphicsItem
{
public:
  ClustersGraphicsItem(const Mesher* mesher)
      : mesher_(mesher)
      , bounds_(-100.0, -100.0, 200.0, 200.0)
      , last_vertex_(DT::Vertex_handle()) {
    setAcceptHoverEvents(true);
    rebuild();
  }

  void setMesher(const Mesher* mesher) {
    mesher_ = mesher;
    rebuild();
    update();
  }

  void setVisible(bool b) {
    if(b) {
      rebuild();
    }
    QGraphicsItem::setVisible(b);
  }

  void rebuild() {
    dt_.clear();
    if(mesher_ == nullptr) {
      return;
    }

    for(auto it : clusters_vertices_range(mesher_->clusters())) {
      dt_.insert(it->point());
    }
  }

  void updateBounds(const QRectF& rect) {
    bounds_ = rect;
    prepareGeometryChange();
    update();
  }

  QRectF boundingRect() const override { return bounds_; }

  void paint(QPainter* painter, const QStyleOptionGraphicsItem*, QWidget*) override {
    if(mesher_ == nullptr) {
      return;
    }

    CGAL::Qt::Converter<K> convert;
    painter->setPen(pen_);
    painter->setBrush(brush_);

    for(auto it : clusters_vertices_range(mesher_->clusters())) {
      QPointF p = convert(it->point());
      painter->drawEllipse(p, 2.0, 2.0);
    }

    if(highlight_segments_.empty()) {
      return;
    }

    QPen reduced_pen(QColor(0, 120, 255), 0);
    QPen regular_pen(QColor(200, 0, 0), 0);
    for(const auto& entry : highlight_segments_) {
      painter->setPen(entry.is_reduced ? reduced_pen : regular_pen);
      painter->drawLine(convert(entry.segment));
    }
  }

  void modelChanged() override {
    // do nothing
  }

  void setPen(const QPen& pen) {
    pen_ = pen;
    update();
  }

  QPen pen() const { return pen_; }

  void setBrush(const QBrush& brush) {
    brush_ = brush;
    update();
  }

  QBrush brush() const { return brush_; }

protected:
  void hoverMoveEvent(QGraphicsSceneHoverEvent* event) override {
    if(mesher_ == nullptr || dt_.number_of_vertices() == 0) {
      clearHighlight();
      update();
      return;
    }

    CGAL::Qt::Converter<K> convert;
    Point_2 p = convert(event->scenePos());
    DT::Vertex_handle nearest = dt_.nearest_vertex(p);
    if(nearest == DT::Vertex_handle()) {
      clearHighlight();
      update();
      return;
    }
    if(nearest == last_vertex_) {
      return;
    }
    last_vertex_ = nearest;

    highlight_segments_.clear();

    typename Mesher::Triangulation::Locate_type lt;
    int i = 0;
    typename Mesher::Triangulation::Face_handle fh = mesher_->triangulation().locate(nearest->point(), lt, i);
    if(lt != Mesher::Triangulation::VERTEX) {
      update();
      return;
    }

    typename Mesher::Triangulation::Vertex_handle v2 = fh->vertex(i);
    int n = mesher_->clusters().number_of_clusters_at_vertex(v2);
    for(int j = 0; j < n; ++j) {
      auto seq = mesher_->clusters().vertices_in_cluster_sequence(v2, j);
      typename Mesher::Clusters::Cluster cluster;
      typename Mesher::Clusters::const_iterator dummy_it;
      mesher_->clusters().get_cluster(v2, *(seq.first), cluster, dummy_it);
      for(auto it = seq.first; it != seq.second; ++it) {
        HighlightSegment entry{Segment_2(v2->point(), (*it)->point()), cluster.is_reduced()};
        highlight_segments_.push_back(entry);
      }
    }
    update();
  }

  void hoverLeaveEvent(QGraphicsSceneHoverEvent*) override {
    clearHighlight();
    update();
  }

private:
  struct HighlightSegment
  {
    Segment_2 segment;
    bool is_reduced;
  };

  void clearHighlight() {
    highlight_segments_.clear();
    last_vertex_ = DT::Vertex_handle();
  }

private:
  using DT = CGAL::Delaunay_triangulation_2<K>;

  const Mesher* mesher_;
  DT dt_;
  QRectF bounds_;
  DT::Vertex_handle last_vertex_;
  std::vector<HighlightSegment> highlight_segments_;
  QPen pen_{QColor(0, 0, 0), 0};
  QBrush brush_{QColor(0, 0, 255, 160)};
};

class LocalSizeRefiner : public CGAL::Qt::GraphicsViewInput
{
public:
  LocalSizeRefiner(QGraphicsScene*, CDT* cdt, Criteria* criteria, Mesher** mesher, QObject* parent)
      : CGAL::Qt::GraphicsViewInput(parent)
      , cdt_(cdt)
      , criteria_(criteria)
      , mesher_(mesher)
      , has_previous_(false) {}

  void reset() { has_previous_ = false; }

protected:
  bool eventFilter(QObject* obj, QEvent* event) override {
    if(event->type() == QEvent::GraphicsSceneMouseMove) {
      QGraphicsSceneMouseEvent* mouseEvent = static_cast<QGraphicsSceneMouseEvent*>(event);
      mouseMoveEvent(mouseEvent);
      return false;
    }
    return QObject::eventFilter(obj, event);
  }

  void mouseMoveEvent(QGraphicsSceneMouseEvent* event) {
    if(cdt_ == nullptr || criteria_ == nullptr || mesher_ == nullptr) {
      return;
    }
    if(*mesher_ == nullptr) {
      return;
    }
    if(cdt_->dimension() != 2) {
      return;
    }

    CGAL::Qt::Converter<K> convert;
    Point_2 current_point = convert(event->scenePos());

    Point_2 previous_point = current_point;
    if(has_previous_) {
      previous_point = previous_point_;
    }

    Face_handle fh = cdt_->locate(current_point);
    if(cdt_->is_infinite(fh)) {
      previous_point_ = current_point;
      has_previous_ = true;
      return;
    }

    std::vector<Face_handle> faces_to_check;
    Segment_2 segment(previous_point, current_point);

    if(has_previous_) {
      auto fc = cdt_->line_walk(current_point, previous_point, fh);
      do {
        faces_to_check.push_back(fc);
        ++fc;
      } while(!cdt_->is_infinite(fc) && cdt_->triangle(fc).has_on_unbounded_side(previous_point));
    } else {
      faces_to_check.push_back(fh);
    }

    criteria_->set_local_size(true);
    criteria_->set_segment(segment);

    std::vector<Face_handle> bad_faces;
    for(const auto& face : faces_to_check) {
      if((face != nullptr) && (!cdt_->is_infinite(face)) && face->is_in_domain()) {
        Criteria::Quality q;
        if(criteria_->is_bad_object().operator()(face, q) != CGAL::Mesh_2::NOT_BAD) {
          bad_faces.push_back(face);
        }
      }
    }

    (*mesher_)->set_criteria(*criteria_, false);
    (*mesher_)->set_bad_faces(bad_faces.begin(), bad_faces.end());
    while((*mesher_)->step_by_step_refine_mesh()) {
    }

    previous_point_ = current_point;
    has_previous_ = true;
  }

private:
  CDT* cdt_;
  Criteria* criteria_;
  Mesher** mesher_;
  Point_2 previous_point_;
  bool has_previous_;
};

class CdtPointInputWithConflictZone : public CGAL::Qt::GraphicsViewInput
{
public:
  CdtPointInputWithConflictZone(QGraphicsScene* scene, CDT* cdt, QObject* parent)
      : CGAL::Qt::GraphicsViewInput(parent)
      , scene_(scene)
      , cdt_(cdt)
      , active_(false) {}

protected:
  bool eventFilter(QObject* obj, QEvent* event) override {
    if(event->type() == QEvent::GraphicsSceneMousePress) {
      QGraphicsSceneMouseEvent* mouseEvent = static_cast<QGraphicsSceneMouseEvent*>(event);
      mousePressEvent(mouseEvent);
      return true;
    }
    if(event->type() == QEvent::GraphicsSceneMouseMove) {
      QGraphicsSceneMouseEvent* mouseEvent = static_cast<QGraphicsSceneMouseEvent*>(event);
      mouseMoveEvent(mouseEvent);
      return false;
    }
    if(event->type() == QEvent::GraphicsSceneMouseRelease) {
      QGraphicsSceneMouseEvent* mouseEvent = static_cast<QGraphicsSceneMouseEvent*>(event);
      mouseReleaseEvent(mouseEvent);
      return true;
    }
    return QObject::eventFilter(obj, event);
  }

  void mousePressEvent(QGraphicsSceneMouseEvent* event) {
    if(event->modifiers() != 0 || event->button() != ::Qt::LeftButton) {
      return;
    }
    if(cdt_ == nullptr || cdt_->number_of_vertices() == 0) {
      return;
    }
    active_ = true;
    last_point_ = convert_(event->scenePos());
    updateConflictZone(event->scenePos());
  }

  void mouseMoveEvent(QGraphicsSceneMouseEvent* event) {
    if(!active_) {
      return;
    }
    last_point_ = convert_(event->scenePos());
    updateConflictZone(event->scenePos());
  }

  void mouseReleaseEvent(QGraphicsSceneMouseEvent*) {
    if(!active_) {
      return;
    }
    clearItems();
    active_ = false;
    Q_EMIT(generate(CGAL::make_object(last_point_)));
  }

private:
  void updateConflictZone(const QPointF& scene_pos) {
    clearItems();
    if(cdt_ == nullptr) {
      return;
    }

    Point_2 p = convert_(scene_pos);
    hint_ = cdt_->locate(p, hint_);

    std::list<CDT::Edge> edges;
    cdt_->find_conflicts(p, edges, hint_);

    QPen pen(QColor(0, 80, 200, 180), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
    for(const auto& edge : edges) {
      Segment_2 s = cdt_->segment(edge);
      QLineF line = convert_(s);
      auto item = new QGraphicsLineItem(line);
      item->setPen(pen);
      item->setZValue(5.0);
      scene_->addItem(item);
      conflict_items_.push_back(item);
    }
  }

  void clearItems() {
    for(auto item : conflict_items_) {
      scene_->removeItem(item);
      delete item;
    }
    conflict_items_.clear();
  }

private:
  QGraphicsScene* scene_;
  CDT* cdt_;
  CGAL::Qt::Converter<K> convert_;
  CDT::Face_handle hint_;
  std::vector<QGraphicsLineItem*> conflict_items_;
  Point_2 last_point_;
  bool active_;
};

class DisplaySettingsWidget : public QWidget
{
  Q_OBJECT

public:
  DisplaySettingsWidget(QWidget* parent)
      : QWidget(parent) {
    QFormLayout* layout = new QFormLayout(this);

    addColorAndWidthRow(layout, "Edges", show_edges_, edge_color_, edge_width_);
    addColorAndWidthRow(layout, "Constraints", show_constraints_, constraint_color_, constraint_width_);
    addColorAndWidthRow(layout, "Vertices", show_vertices_, vertex_color_, vertex_width_);
    addColorOnlyRow(layout, "Faces in domain", show_faces_, faces_color_);
    addColorAndWidthRow(layout, "Seeds", show_seeds_, seeds_color_, seeds_width_);
    addColorOnlyRow(layout, "Bad faces", show_bad_faces_, bad_faces_color_);
    addColorAndWidthRow(layout, "Encroached edges", show_encroached_edges_, encroached_edges_color_, encroached_edges_width_);
    addColorAndWidthRow(layout, "Clusters", show_clusters_, clusters_color_, clusters_width_);

    setLayout(layout);
  }



  QColor edgeColor() const { return buttonColor(edge_color_); }
  QColor constraintColor() const { return buttonColor(constraint_color_); }
  QColor vertexColor() const { return buttonColor(vertex_color_); }
  QColor facesColor() const { return buttonColor(faces_color_); }
  QColor seedsColor() const { return buttonColor(seeds_color_); }
  QColor badFacesColor() const { return buttonColor(bad_faces_color_); }
  QColor encroachedEdgesColor() const { return buttonColor(encroached_edges_color_); }
  QColor clustersColor() const { return buttonColor(clusters_color_); }

  int edgeWidth() const { return edge_width_->value(); }
  int constraintWidth() const { return constraint_width_->value(); }
  int vertexWidth() const { return vertex_width_->value(); }
  int seedsWidth() const { return seeds_width_->value(); }
  int encroachedEdgesWidth() const { return encroached_edges_width_->value(); }
  int clustersWidth() const { return clusters_width_->value(); }



  bool showEdges() const { return show_edges_->isChecked(); }
  bool showConstraints() const { return show_constraints_->isChecked(); }
  bool showVertices() const { return show_vertices_->isChecked(); }
  bool showFaces() const { return show_faces_->isChecked(); }
  bool showSeeds() const { return show_seeds_->isChecked(); }
  bool showBadFaces() const { return show_bad_faces_->isChecked(); }
  bool showEncroachedEdges() const { return show_encroached_edges_->isChecked(); }
  bool showClusters() const { return show_clusters_->isChecked(); }

  void saveSettings() const {
    QSettings settings("GeometryFactory", "Mesh_2_demo");
    settings.beginGroup("DisplaySettings");

    settings.setValue("edgeColor", edgeColor());
    settings.setValue("edgeWidth", edgeWidth());
    settings.setValue("constraintColor", constraintColor());
    settings.setValue("constraintWidth", constraintWidth());
    settings.setValue("vertexColor", vertexColor());
    settings.setValue("vertexWidth", vertexWidth());
    settings.setValue("facesColor", facesColor());
    settings.setValue("seedsColor", seedsColor());
    settings.setValue("seedsWidth", seedsWidth());
    settings.setValue("badFacesColor", badFacesColor());
    settings.setValue("encroachedEdgesColor", encroachedEdgesColor());
    settings.setValue("encroachedEdgesWidth", encroachedEdgesWidth());
    settings.setValue("clustersColor", clustersColor());
    settings.setValue("clustersWidth", clustersWidth());

    settings.endGroup();
  }

  void restoreSettings() {
    QSettings settings("GeometryFactory", "Mesh_2_demo");
    settings.beginGroup("DisplaySettings");

    setButtonColor(edge_color_, settings.value("edgeColor", QColor(Qt::blue)).value<QColor>());
    edge_width_->setValue(settings.value("edgeWidth", 0).toInt());
    setButtonColor(constraint_color_, settings.value("constraintColor", QColor(Qt::red)).value<QColor>());
    constraint_width_->setValue(settings.value("constraintWidth", 0).toInt());
    setButtonColor(vertex_color_, settings.value("vertexColor", QColor(Qt::black)).value<QColor>());
    vertex_width_->setValue(settings.value("vertexWidth", 2).toInt());
    setButtonColor(faces_color_, settings.value("facesColor", QColor(Qt::green)).value<QColor>());
    setButtonColor(seeds_color_, settings.value("seedsColor", QColor(Qt::darkBlue)).value<QColor>());
    seeds_width_->setValue(settings.value("seedsWidth", 0).toInt());
    setButtonColor(bad_faces_color_, settings.value("badFacesColor", QColor(0, 160, 0, 128)).value<QColor>());
    setButtonColor(encroached_edges_color_, settings.value("encroachedEdgesColor", QColor(Qt::red)).value<QColor>());
    encroached_edges_width_->setValue(settings.value("encroachedEdgesWidth", 0).toInt());
    setButtonColor(clusters_color_, settings.value("clustersColor", QColor(Qt::black)).value<QColor>());
    clusters_width_->setValue(settings.value("clustersWidth", 0).toInt());

    settings.endGroup();
  }

Q_SIGNALS:
  void displayChanged();

private:
  void addColorAndWidthRow(QFormLayout* layout, const QString& label, QCheckBox*& checkbox, QPushButton*& colorButton, QSpinBox*& widthSpinBox) {
    QHBoxLayout* rowLayout = new QHBoxLayout();

    checkbox = new QCheckBox(this);
    checkbox->setChecked(true);
    connect(checkbox, &QCheckBox::toggled, this, &DisplaySettingsWidget::displayChanged);

    colorButton = new QPushButton(this);
    colorButton->setAutoFillBackground(true);
    colorButton->setMinimumWidth(80);
    connect(colorButton, &QPushButton::clicked, this, &DisplaySettingsWidget::selectColor);

    widthSpinBox = new QSpinBox(this);
    widthSpinBox->setRange(0, 100);
    widthSpinBox->setValue(0);
    widthSpinBox->setPrefix("Width: ");
    connect(widthSpinBox, QOverload<int>::of(&QSpinBox::valueChanged), this, &DisplaySettingsWidget::displayChanged);

    rowLayout->addWidget(checkbox);
    rowLayout->addWidget(colorButton);
    rowLayout->addWidget(widthSpinBox);
    rowLayout->addStretch();

    layout->addRow(label, rowLayout);
  }

  void addColorOnlyRow(QFormLayout* layout, const QString& label, QCheckBox*& checkbox, QPushButton*& colorButton) {
    QHBoxLayout* rowLayout = new QHBoxLayout();

    checkbox = new QCheckBox(this);
    checkbox->setChecked(true);
    connect(checkbox, &QCheckBox::toggled, this, &DisplaySettingsWidget::displayChanged);

    colorButton = new QPushButton(this);
    colorButton->setAutoFillBackground(true);
    colorButton->setMinimumWidth(80);
    connect(colorButton, &QPushButton::clicked, this, &DisplaySettingsWidget::selectColor);

    rowLayout->addWidget(checkbox);
    rowLayout->addWidget(colorButton);
    rowLayout->addStretch();

    layout->addRow(label, rowLayout);
  }

  void setButtonColor(QPushButton* button, const QColor& color) {
    std::cerr << "Setting button color to " << color.name().toStdString() << std::endl;
    QPalette pal = button->palette();
    pal.setColor(QPalette::Button, color);
    button->setPalette(pal);
    button->setText(color.name());
    button->update();
  }

  QColor buttonColor(QPushButton* button) const { return button->palette().color(QPalette::Button); }

private Q_SLOTS:
  void selectColor() {
    QPushButton* button = qobject_cast<QPushButton*>(sender());
    if(button == nullptr) {
      return;
    }
    QColor color = QColorDialog::getColor(buttonColor(button), this);
    if(color.isValid()) {
      setButtonColor(button, color);
      Q_EMIT displayChanged();
    }
  }

private:
  QCheckBox* show_edges_;
  QPushButton* edge_color_;
  QSpinBox* edge_width_;
  QCheckBox* show_constraints_;
  QPushButton* constraint_color_;
  QSpinBox* constraint_width_;
  QCheckBox* show_vertices_;
  QPushButton* vertex_color_;
  QSpinBox* vertex_width_;
  QCheckBox* show_faces_;
  QPushButton* faces_color_;
  QCheckBox* show_seeds_;
  QPushButton* seeds_color_;
  QSpinBox* seeds_width_;
  QCheckBox* show_bad_faces_;
  QPushButton* bad_faces_color_;
  QCheckBox* show_encroached_edges_;
  QPushButton* encroached_edges_color_;
  QSpinBox* encroached_edges_width_;
  QCheckBox* show_clusters_;
  QPushButton* clusters_color_;
  QSpinBox* clusters_width_;
};

class MainWindow : public CGAL::Qt::DemosMainWindow
{
  Q_OBJECT

public:
  MainWindow();
  void openTriangulation(const QString& filename);

private Q_SLOTS:
  void onActionInsertPolylineToggled(bool checked);
  void onActionInsertPointsToggled(bool checked);
  void onActionInsertSeedsToggled(bool checked);
  void onActionShowCircumcircleToggled(bool checked);

  void onActionInsertBoundingBox();
  void onActionRefineMesh();
  void onActionConformMesh();
  void onActionRefineStep();
  void onActionAutoStepToggled(bool checked);
  void onActionClear();
  void onActionRecenter();
  void onActionOpen();
  void onActionSave();


  void processInput(CGAL::Object o);
  void updateStepLength(int value);
  void updateTimerInterval(int value);
  void onUnderMouseToggled(bool checked);
  void applyDisplaySettings();
  void updateCriteria();

private:
  void setupUi();
  void setupScene();
  void setupActions();
  void setupDock();
  void updateCounts();
  void updateMesherState();
  void updateDomain();
  void updateBounds();
  void restoreVisibilitySettings();
  void saveVisibilitySettings();
  void saveTriangulation(const QString& filename);
  void loadPolyConstraints(const QString& filename);
  void loadEdgConstraints(const QString& filename);
  void loadDataConstraints(const QString& filename);
  void loadObjConstraints(const QString& filename);

  Mesher* createMesher();

private:
  CDT cdt_;
  Criteria criteria_;
  Mesher* mesher_ = nullptr;
  std::list<Point_2> seeds_;

  QGraphicsScene scene_;
  QGraphicsView* view_ = nullptr;

  CGAL::Qt::DelaunayMeshTriangulationGraphicsItem<CDT>* dgi_ = nullptr;
  CGAL::Qt::GraphicsViewPolylineInput<K>* polyline_input_ = nullptr;
  CdtPointInputWithConflictZone* point_input_ = nullptr;
  CGAL::Qt::DelaunayMeshInsertSeeds<CDT>* seed_input_ = nullptr;
  CGAL::Qt::TriangulationCircumcircle<CDT>* circumcircle_ = nullptr;
  LocalSizeRefiner* local_refiner_ = nullptr;

  BadFacesGraphicsItem* bad_faces_item_ = nullptr;
  EncroachedEdgesGraphicsItem* encroached_edges_item_ = nullptr;
  ClustersGraphicsItem* clusters_item_ = nullptr;

  QAction* action_show_vertices_ = nullptr;
  QAction* action_show_edges_ = nullptr;
  QAction* action_show_constraints_ = nullptr;
  QAction* action_show_faces_ = nullptr;
  QAction* action_show_seeds_ = nullptr;
  QAction* action_show_bad_faces_ = nullptr;
  QAction* action_show_encroached_edges_ = nullptr;
  QAction* action_show_clusters_ = nullptr;
  QAction* action_auto_step_ = nullptr;

  QLabel* nb_points_label_ = nullptr;
  QLabel* nb_clusters_label_ = nullptr;
  QLabel* init_status_label_ = nullptr;
  QDoubleSpinBox* angle_bound_box_ = nullptr;
  QDoubleSpinBox* size_bound_box_ = nullptr;
  QCheckBox* under_mouse_box_ = nullptr;
  QSpinBox* step_length_box_ = nullptr;
  QSpinBox* timer_interval_box_ = nullptr;

  QTimer* timer_ = nullptr;
  int step_length_ = 1;
  int timer_interval_ = 1000;

  DisplaySettingsWidget* display_settings_ = nullptr;
  QDockWidget* display_dock_ = nullptr;
};

MainWindow::MainWindow() {
  setupUi();
  setupScene();
  setupActions();
  setupDock();

  // Restore display settings after dock is set up
  display_settings_->restoreSettings();
  connect(display_settings_, &DisplaySettingsWidget::displayChanged, this, &MainWindow::applyDisplaySettings);
  applyDisplaySettings();
  restoreVisibilitySettings();

  onActionInsertPolylineToggled(true);
  dgi_->setVisibleSeeds(action_show_seeds_->isChecked(), seeds_.begin(), seeds_.end());
  updateCriteria();
  updateCounts();
}

void MainWindow::setupUi() {
  view_ = new QGraphicsView(this);
  this->view = view_;
  setCentralWidget(view_);

  statusBar();
  setupStatusBar();
  setupOptionsMenu();
  addAboutDemo(":/cgal/help/about_Constrained_Delaunay_triangulation_2.html");
  addAboutCGAL();
}

void MainWindow::setupScene() {
  scene_.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene_.setSceneRect(-100.0, -100.0, 200.0, 200.0);
  view_->setScene(&scene_);
  view_->setMouseTracking(true);
  view_->scale(1.0, -1.0);
  addNavigation(view_);

  bad_faces_item_ = new BadFacesGraphicsItem(&cdt_, mesher_);
  bad_faces_item_->setVisible(false);
  scene_.addItem(bad_faces_item_);

  encroached_edges_item_ = new EncroachedEdgesGraphicsItem(&cdt_, mesher_);
  encroached_edges_item_->setVisible(false);
  scene_.addItem(encroached_edges_item_);

  clusters_item_ = new ClustersGraphicsItem(mesher_);
  clusters_item_->setVisible(false);
  scene_.addItem(clusters_item_);

  dgi_ = new CGAL::Qt::DelaunayMeshTriangulationGraphicsItem<CDT>(&cdt_);
  QColor faces_color(Qt::green);
  faces_color.setAlpha(150);
  dgi_->setFacesInDomainBrush(faces_color);
  dgi_->setVerticesPen(pen(Qt::black, 2));
  dgi_->setConstraintsPen(pen(Qt::red));
  dgi_->setEdgesPen(pen(Qt::blue));
  dgi_->setSeedsPen(pen(Qt::darkBlue));
  dgi_->setZValue(-1);
  scene_.addItem(dgi_);

  polyline_input_ = new CGAL::Qt::GraphicsViewPolylineInput<K>(this, &scene_, 0, true);
  connect(polyline_input_, &CGAL::Qt::GraphicsViewInput::generate, this, &MainWindow::processInput);

  point_input_ = new CdtPointInputWithConflictZone(&scene_, &cdt_, this);
  connect(point_input_, &CGAL::Qt::GraphicsViewInput::generate, this, &MainWindow::processInput);

  seed_input_ = new CGAL::Qt::DelaunayMeshInsertSeeds<CDT>(&scene_, &cdt_, this);
  connect(seed_input_, &CGAL::Qt::GraphicsViewInput::generate, this, &MainWindow::processInput);

  circumcircle_ = new CGAL::Qt::TriangulationCircumcircle<CDT>(&scene_, &cdt_, this);
  circumcircle_->setPen(pen(Qt::darkGray));

  local_refiner_ = new LocalSizeRefiner(&scene_, &cdt_, &criteria_, &mesher_, this);
}

void MainWindow::setupActions() {
  QMenu* file_menu = menuBar()->addMenu("&File");
  QMenu* mesh_menu = menuBar()->addMenu("&Mesh");
  QMenu* view_menu = menuBar()->addMenu("&View");
  QMenu* input_menu = menuBar()->addMenu("&Input");
  QMenu* options_menu = menuBar()->addMenu("&Options");

  QAction* action_open = file_menu->addAction("Open constraints...");
  QAction* action_save = file_menu->addAction("Save constraints...");
  QAction* action_clear = file_menu->addAction("Clear");
  QAction* action_recenter = file_menu->addAction("Recenter");
  file_menu->addSeparator();
  QAction* action_quit = file_menu->addAction("Quit");

  QAction* action_insert_bounding_box = mesh_menu->addAction("Insert bounding box");
  QAction* action_refine_mesh = mesh_menu->addAction("Refine mesh");
  QAction* action_conform_mesh = mesh_menu->addAction("Conform mesh");
  QAction* action_refine_step = mesh_menu->addAction("Refine one step");
  action_auto_step_ = mesh_menu->addAction("Auto step");
  action_auto_step_->setCheckable(true);
  QAction* action_auto_step = action_auto_step_;

  action_show_vertices_ = view_menu->addAction("Show vertices");
  action_show_edges_ = view_menu->addAction("Show edges");
  action_show_constraints_ = view_menu->addAction("Show constrained edges");
  action_show_faces_ = view_menu->addAction("Show faces in domain");
  action_show_seeds_ = view_menu->addAction("Show seeds");
  action_show_bad_faces_ = view_menu->addAction("Show bad faces");
  action_show_encroached_edges_ = view_menu->addAction("Show encroached edges");
  action_show_clusters_ = view_menu->addAction("Show clusters");

  action_show_vertices_->setCheckable(true);
  action_show_edges_->setCheckable(true);
  action_show_constraints_->setCheckable(true);
  action_show_faces_->setCheckable(true);
  action_show_seeds_->setCheckable(true);
  action_show_bad_faces_->setCheckable(true);
  action_show_encroached_edges_->setCheckable(true);
  action_show_clusters_->setCheckable(true);

  action_show_vertices_->setChecked(true);
  action_show_edges_->setChecked(true);
  action_show_constraints_->setChecked(true);
  action_show_faces_->setChecked(true);
  action_show_seeds_->setChecked(true);

  QAction* action_insert_polyline = input_menu->addAction("Insert polyline");
  QAction* action_insert_points = input_menu->addAction("Insert points");
  QAction* action_insert_seeds = input_menu->addAction("Insert seeds");
  QAction* action_show_circumcircle = input_menu->addAction("Show circumcircle");

  action_insert_polyline->setCheckable(true);
  action_insert_points->setCheckable(true);
  action_insert_seeds->setCheckable(true);
  action_show_circumcircle->setCheckable(true);

  QActionGroup* input_group = new QActionGroup(this);
  input_group->addAction(action_insert_polyline);
  input_group->addAction(action_insert_points);
  input_group->addAction(action_insert_seeds);
  action_insert_polyline->setChecked(true);

  QAction* action_display_settings = options_menu->addAction("Display settings...");

  QToolBar* toolbar = addToolBar("Mesh actions");
  toolbar->addAction(action_open);
  toolbar->addAction(action_save);
  toolbar->addAction(action_refine_mesh);
  toolbar->addAction(action_conform_mesh);
  toolbar->addAction(action_refine_step);
  toolbar->addAction(action_auto_step);

  QToolBar* input_toolbar = addToolBar("Input");
  input_toolbar->addAction(action_insert_polyline);
  input_toolbar->addAction(action_insert_points);
  input_toolbar->addAction(action_insert_seeds);
  input_toolbar->addAction(action_show_circumcircle);

  connect(action_open, &QAction::triggered, this, &MainWindow::onActionOpen);
  connect(action_save, &QAction::triggered, this, &MainWindow::onActionSave);
  connect(action_clear, &QAction::triggered, this, &MainWindow::onActionClear);
  connect(action_recenter, &QAction::triggered, this, &MainWindow::onActionRecenter);
  connect(action_quit, &QAction::triggered, this, &MainWindow::close);

  connect(action_insert_bounding_box, &QAction::triggered, this, &MainWindow::onActionInsertBoundingBox);
  connect(action_refine_mesh, &QAction::triggered, this, &MainWindow::onActionRefineMesh);
  connect(action_conform_mesh, &QAction::triggered, this, &MainWindow::onActionConformMesh);
  connect(action_refine_step, &QAction::triggered, this, &MainWindow::onActionRefineStep);
  connect(action_auto_step_, &QAction::toggled, this, &MainWindow::onActionAutoStepToggled);

  using DtGraphicsItem = CGAL::Qt::DelaunayMeshTriangulationGraphicsItem<CDT>;

  connect(action_show_vertices_, &QAction::toggled, dgi_, &DtGraphicsItem::setVisibleVertices);
  connect(action_show_edges_, &QAction::toggled, dgi_, &DtGraphicsItem::setVisibleEdges);
  connect(action_show_constraints_, &QAction::toggled, dgi_, &DtGraphicsItem::setVisibleConstraints);
  connect(action_show_faces_, &QAction::toggled, dgi_, &DtGraphicsItem::setVisibleFacesInDomain);
  connect(action_show_seeds_, &QAction::toggled,
          [this](bool checked) { dgi_->setVisibleSeeds(checked, seeds_.begin(), seeds_.end()); });

  using CGAL_graphicsItem = CGAL::Qt::GraphicsItem;
  connect(action_show_bad_faces_, &QAction::toggled, bad_faces_item_, &CGAL_graphicsItem::setVisible);
  connect(action_show_encroached_edges_, &QAction::toggled, encroached_edges_item_, &CGAL_graphicsItem::setVisible);
  connect(action_show_clusters_, &QAction::toggled, clusters_item_, &CGAL_graphicsItem::setVisible);

  connect(action_show_vertices_, &QAction::toggled, this, [this](bool) { saveVisibilitySettings(); });
  connect(action_show_edges_, &QAction::toggled, this, [this](bool) { saveVisibilitySettings(); });
  connect(action_show_constraints_, &QAction::toggled, this, [this](bool) { saveVisibilitySettings(); });
  connect(action_show_faces_, &QAction::toggled, this, [this](bool) { saveVisibilitySettings(); });
  connect(action_show_seeds_, &QAction::toggled, this, [this](bool) { saveVisibilitySettings(); });
  connect(action_show_bad_faces_, &QAction::toggled, this, [this](bool) { saveVisibilitySettings(); });
  connect(action_show_encroached_edges_, &QAction::toggled, this, [this](bool) { saveVisibilitySettings(); });
  connect(action_show_clusters_, &QAction::toggled, this, [this](bool) { saveVisibilitySettings(); });

  connect(action_insert_polyline, &QAction::toggled, this, &MainWindow::onActionInsertPolylineToggled);
  connect(action_insert_points, &QAction::toggled, this, &MainWindow::onActionInsertPointsToggled);
  connect(action_insert_seeds, &QAction::toggled, this, &MainWindow::onActionInsertSeedsToggled);
  connect(action_show_circumcircle, &QAction::toggled, this, &MainWindow::onActionShowCircumcircleToggled);

  connect(action_display_settings, &QAction::triggered, [this]() {
    display_dock_->setVisible(!display_dock_->isVisible());
  });
}

void MainWindow::setupDock() {
  QDockWidget* dock = new QDockWidget("Mesh settings", this);
  QWidget* dock_widget = new QWidget(dock);
  QVBoxLayout* layout = new QVBoxLayout(dock_widget);

  QGroupBox* info_box = new QGroupBox("Info", dock_widget);
  QFormLayout* info_layout = new QFormLayout(info_box);
  nb_points_label_ = new QLabel("0", info_box);
  nb_clusters_label_ = new QLabel("0", info_box);
  init_status_label_ = new QLabel("no", info_box);
  info_layout->addRow("Number of points", nb_points_label_);
  info_layout->addRow("Number of clusters", nb_clusters_label_);
  info_layout->addRow("Initialized", init_status_label_);

  QGroupBox* criteria_box = new QGroupBox("Criteria", dock_widget);
  QFormLayout* criteria_layout = new QFormLayout(criteria_box);
  angle_bound_box_ = new QDoubleSpinBox(criteria_box);
  angle_bound_box_->setDecimals(4);
  angle_bound_box_->setRange(0.0, 1000.0);
  angle_bound_box_->setValue(0.125);
  size_bound_box_ = new QDoubleSpinBox(criteria_box);
  size_bound_box_->setDecimals(6);
  size_bound_box_->setRange(0.0, std::numeric_limits<double>::max());
  size_bound_box_->setValue(0.0);
  criteria_layout->addRow("Angle bound", angle_bound_box_);
  criteria_layout->addRow("Size bound", size_bound_box_);

  under_mouse_box_ = new QCheckBox("Under mouse only", criteria_box);
  criteria_layout->addRow(under_mouse_box_);

  QGroupBox* step_box = new QGroupBox("Step controls", dock_widget);
  QFormLayout* step_layout = new QFormLayout(step_box);
  step_length_box_ = new QSpinBox(step_box);
  step_length_box_->setRange(1, std::numeric_limits<int>::max());
  step_length_box_->setValue(step_length_);
  timer_interval_box_ = new QSpinBox(step_box);
  timer_interval_box_->setRange(1, std::numeric_limits<int>::max());
  timer_interval_box_->setValue(timer_interval_);
  timer_interval_box_->setSuffix(" ms");
  step_layout->addRow("Step length", step_length_box_);
  step_layout->addRow("Timer interval", timer_interval_box_);

  layout->addWidget(info_box);
  layout->addWidget(criteria_box);
  layout->addWidget(step_box);
  layout->addStretch(1);

  dock_widget->setLayout(layout);
  dock->setWidget(dock_widget);
  addDockWidget(Qt::RightDockWidgetArea, dock);

  connect(angle_bound_box_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::updateCriteria);
  connect(size_bound_box_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::updateCriteria);
  connect(under_mouse_box_, &QCheckBox::toggled, this, &MainWindow::onUnderMouseToggled);
  connect(step_length_box_, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::updateStepLength);
  connect(timer_interval_box_, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::updateTimerInterval);

  // Create display settings dock widget
  display_dock_ = new QDockWidget("Display settings", this);
  display_settings_ = new DisplaySettingsWidget(display_dock_);
  display_dock_->setWidget(display_settings_);
  addDockWidget(Qt::RightDockWidgetArea, display_dock_);

  timer_ = new QTimer(this);
  connect(timer_, &QTimer::timeout, this, &MainWindow::onActionRefineStep);
}

void MainWindow::processInput(CGAL::Object o) {
  std::list<Point_2> points;
  if(CGAL::assign(points, o)) {
    if(points.size() == 1) {
      cdt_.insert(points.front());
    } else {
      Point_2 p, q;
      CDT::Vertex_handle vh, wh;
      auto it = points.begin();
      vh = cdt_.insert(*it);
      p = *it;
      ++it;
      for(; it != points.end(); ++it) {
        q = *it;
        if(p != q) {
          wh = cdt_.insert(*it);
          cdt_.insert_constraint(vh, wh);
          vh = wh;
          p = q;
        }
      }
    }
  } else {
    Point_2 p;
    if(CGAL::assign(p, o)) {
      seeds_.push_back(p);
      if(action_show_seeds_->isChecked()) {
        dgi_->setVisibleSeeds(true, seeds_.begin(), seeds_.end());
      }
    }
  }

  updateDomain();
  updateCounts();
  updateBounds();
  dgi_->modelChanged();
}

void MainWindow::updateDomain() {
  Mesher::mark_facets(cdt_, seeds_.begin(), seeds_.end());
}

void MainWindow::updateCounts() {
  nb_points_label_->setText(QString::number(static_cast<int>(cdt_.number_of_vertices())));
  if(mesher_ != nullptr) {
    nb_clusters_label_->setText(QString::number(mesher_->clusters().size()));
    init_status_label_->setText("yes");
  } else {
    nb_clusters_label_->setText("0");
    init_status_label_->setText("no");
  }
}

void MainWindow::updateMesherState() {
  bad_faces_item_->setMesher(mesher_);
  encroached_edges_item_->setMesher(mesher_);
  clusters_item_->setMesher(mesher_);
}

void MainWindow::updateBounds() {
  QRectF rect = dgi_->boundingRect();
  if(!rect.isValid()) {
    rect = QRectF(-100.0, -100.0, 200.0, 200.0);
  }
  bad_faces_item_->updateBounds(rect);
  encroached_edges_item_->updateBounds(rect);
  clusters_item_->updateBounds(rect);
}

void MainWindow::restoreVisibilitySettings() {
  QSettings settings("GeometryFactory", "Mesh_2_demo");
  std::cerr << "Restoring visibility settings from file \"" << settings.fileName().toStdString() << "\"...\n";
  settings.beginGroup("DisplayVisibility");

  action_show_vertices_->setChecked(settings.value("showVertices", true).toBool());
  action_show_edges_->setChecked(settings.value("showEdges", true).toBool());
  action_show_constraints_->setChecked(settings.value("showConstraints", true).toBool());
  action_show_faces_->setChecked(settings.value("showFaces", true).toBool());
  action_show_seeds_->setChecked(settings.value("showSeeds", true).toBool());
  action_show_bad_faces_->setChecked(settings.value("showBadFaces", false).toBool());
  action_show_encroached_edges_->setChecked(settings.value("showEncroachedEdges", false).toBool());
  action_show_clusters_->setChecked(settings.value("showClusters", false).toBool());

  settings.endGroup();
}

void MainWindow::saveVisibilitySettings() {
  QSettings settings("GeometryFactory", "Mesh_2_demo");
  std::cerr << "Saving visibility settings to file \"" << settings.fileName().toStdString() << "\"...\n";
  settings.beginGroup("DisplayVisibility");

  settings.setValue("showVertices", action_show_vertices_->isChecked());
  settings.setValue("showEdges", action_show_edges_->isChecked());
  settings.setValue("showConstraints", action_show_constraints_->isChecked());
  settings.setValue("showFaces", action_show_faces_->isChecked());
  settings.setValue("showSeeds", action_show_seeds_->isChecked());
  settings.setValue("showBadFaces", action_show_bad_faces_->isChecked());
  settings.setValue("showEncroachedEdges", action_show_encroached_edges_->isChecked());
  settings.setValue("showClusters", action_show_clusters_->isChecked());

  settings.endGroup();
}

Mesher* MainWindow::createMesher() {
  Mesher* mesher = new Mesher(cdt_, criteria_);
  mesher->set_seeds(seeds_.begin(), seeds_.end());
  return mesher;
}

void MainWindow::onActionInsertPolylineToggled(bool checked) {
  if(checked) {
    scene_.installEventFilter(polyline_input_);
  } else {
    scene_.removeEventFilter(polyline_input_);
  }
}

void MainWindow::onActionInsertPointsToggled(bool checked) {
  if(checked) {
    scene_.installEventFilter(point_input_);
  } else {
    scene_.removeEventFilter(point_input_);
  }
}

void MainWindow::onActionInsertSeedsToggled(bool checked) {
  if(checked) {
    scene_.installEventFilter(seed_input_);
  } else {
    scene_.removeEventFilter(seed_input_);
  }
}

void MainWindow::onActionShowCircumcircleToggled(bool checked) {
  if(checked) {
    scene_.installEventFilter(circumcircle_);
    circumcircle_->show();
  } else {
    scene_.removeEventFilter(circumcircle_);
    circumcircle_->hide();
  }
}

void MainWindow::onActionInsertBoundingBox() {
  if(cdt_.dimension() < 1) {
    return;
  }

  FT xmin, xmax, ymin, ymax;
  auto vi = cdt_.finite_vertices_begin();
  xmin = xmax = vi->point().x();
  ymin = ymax = vi->point().y();
  ++vi;
  while(vi != cdt_.finite_vertices_end()) {
    xmin = (std::min)(xmin, vi->point().x());
    xmax = (std::max)(xmax, vi->point().x());
    ymin = (std::min)(ymin, vi->point().y());
    ymax = (std::max)(ymax, vi->point().y());
    ++vi;
  }

  FT xcenter = (xmin + xmax) / 2.;
  FT ycenter = (ymin + ymax) / 2.;
  FT xspan = (xmax - xmin) / 2.;
  FT yspan = (ymax - ymin) / 2.;

  Point_2 bb1(xcenter - 1.5 * xspan, ycenter - 1.5 * yspan);
  Point_2 bb2(xcenter + 1.5 * xspan, ycenter - 1.5 * yspan);
  Point_2 bb3(xcenter + 1.5 * xspan, ycenter + 1.5 * yspan);
  Point_2 bb4(xcenter - 1.5 * xspan, ycenter + 1.5 * yspan);

  cdt_.insert(bb1);
  cdt_.insert(bb2);
  cdt_.insert(bb3);
  cdt_.insert(bb4);
  cdt_.insert_constraint(bb1, bb2);
  cdt_.insert_constraint(bb2, bb3);
  cdt_.insert_constraint(bb3, bb4);
  cdt_.insert_constraint(bb4, bb1);

  updateDomain();
  updateCounts();
  dgi_->modelChanged();
}

void MainWindow::onActionRefineMesh() {
  criteria_.set_bound(angle_bound_box_->value());
  criteria_.set_size_bound(size_bound_box_->value());
  criteria_.set_local_size(false);

  if(mesher_ == nullptr) {
    mesher_ = createMesher();
  }

  mesher_->set_criteria(criteria_);
  mesher_->refine_mesh();
  updateMesherState();
  updateDomain();
  updateCounts();
  dgi_->modelChanged();
}

void MainWindow::onActionConformMesh() {
  CGAL::make_conforming_Gabriel_2(cdt_);
  updateDomain();
  if(mesher_ != nullptr) {
    delete mesher_;
    mesher_ = nullptr;
  }
  updateMesherState();
  updateCounts();
  dgi_->modelChanged();
}

void MainWindow::onActionRefineStep() {
  criteria_.set_bound(angle_bound_box_->value());
  criteria_.set_size_bound(size_bound_box_->value());

  if(mesher_ == nullptr) {
    mesher_ = createMesher();
    mesher_->init();
    updateMesherState();
  } else {
    int counter = step_length_;
    while(counter > 0) {
      --counter;
      if(!mesher_->try_one_step_refine_mesh()) {
        action_auto_step_->setChecked(false);
        break;
      }
    }
  }

  updateDomain();
  updateCounts();
  dgi_->modelChanged();
}

void MainWindow::onActionAutoStepToggled(bool checked) {
  if(checked) {
    timer_->start(timer_interval_);
  } else {
    timer_->stop();
  }
}

void MainWindow::onActionClear() {
  cdt_.clear();
  seeds_.clear();
  if(mesher_ != nullptr) {
    delete mesher_;
    mesher_ = nullptr;
  }
  local_refiner_->reset();
  updateMesherState();
  updateCounts();
  dgi_->modelChanged();
}

void MainWindow::onActionRecenter() {
  QRectF rect = dgi_->boundingRect();
  if(rect.isValid()) {
    view_->setSceneRect(rect);
    view_->fitInView(rect, Qt::KeepAspectRatio);
  }
}

void MainWindow::onActionOpen() {
  QString filename = QFileDialog::getOpenFileName(this, "Open constraint file", ".",
                                                  "Edge files (*.edg);;"
                                                  "Poly files (*.poly);;"
                                                  "Data files (*.data);;"
                                                  "Obj files (*.obj);;"
                                                  "All files (*)");
  if(filename.isEmpty()) {
    return;
  }
  openTriangulation(filename);
}

void MainWindow::onActionSave() {
  QString filename = QFileDialog::getSaveFileName(this, "Save constraints", "filename.edg",
                                                  "Edge files (*.edg);;"
                                                  "Poly files (*.poly);;"
                                                  "All files (*)");
  if(filename.isEmpty()) {
    return;
  }
  saveTriangulation(filename);
}

void MainWindow::openTriangulation(const QString& filename) {
  if(filename.endsWith(".poly")) {
    loadPolyConstraints(filename);
  } else if(filename.endsWith(".data")) {
    loadDataConstraints(filename);
  } else if(filename.endsWith(".obj")) {
    loadObjConstraints(filename);
  } else {
    loadEdgConstraints(filename);
  }

  updateDomain();
  updateCounts();
  dgi_->modelChanged();
  onActionRecenter();
}

void MainWindow::saveTriangulation(const QString& filename) {
  std::ofstream of(qPrintable(filename));
  if(!of) {
    return;
  }

  if(filename.endsWith(".poly")) {
    CGAL::IO::write_triangle_poly_file(cdt_, of, seeds_.begin(), seeds_.end());
  } else {
    write_constraints(cdt_, of);
  }
}

void MainWindow::loadPolyConstraints(const QString& filename) {
  std::ifstream ifs(qPrintable(filename));
  if(!ifs) {
    return;
  }
  cdt_.clear();
  seeds_.clear();
  CGAL::IO::read_triangle_poly_file(cdt_, ifs, std::back_inserter(seeds_));
}

void MainWindow::loadEdgConstraints(const QString& filename) {
  std::ifstream ifs(qPrintable(filename));
  if(!ifs) {
    return;
  }
  cdt_.clear();
  seeds_.clear();
  read_constraints(cdt_, ifs);
}

void MainWindow::loadObjConstraints(const QString& filename) {
  std::ifstream ifs(qPrintable(filename));
  if(!ifs) {
    return;
  }
  cdt_.clear();
  seeds_.clear();

  std::vector<Point_2> points;
  std::vector<std::vector<std::size_t>> polylines;
  std::vector<std::vector<std::size_t>> unused_polygons;

  if(!CGAL::IO::internal::read_OBJ(ifs, points, polylines, unused_polygons)) {
    QMessageBox::warning(this, "Error", "Failed to read OBJ file");
    return;
  }

  // Insert constraints from polylines
  for(const std::vector<std::size_t>& polyline : polylines) {
    for(std::size_t i = 1; i < polyline.size(); ++i) {
      cdt_.insert_constraint(points[polyline[i - 1]], points[polyline[i]]);
    }
  }
}

void MainWindow::loadDataConstraints(const QString& filename) {
  std::ifstream ins(qPrintable(filename));
  if(!ins) {
    return;
  }

  int nx, ny, niso, use_threshold;
  float threshold;
  ins >> nx >> ny >> niso >> use_threshold >> threshold;
  for(int c = 0; c < niso; ++c) {
    float f;
    ins >> f;
  }

  std::vector<Point_2> points(nx * ny);
  double xmin, xmax, ymin, ymax;
  ins >> xmin >> xmax >> ymin >> ymax;

  double dx = (xmax - xmin) / (nx - 1);
  double dy = (ymax - ymin) / (ny - 1);

  std::size_t k2 = 0;
  for(int i2 = 0; i2 < nx; ++i2) {
    for(int j = 0; j < ny; ++j) {
      points[k2] = Point_2(xmin + i2 * dx, ymin + j * dy);
      ++k2;
    }
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::shuffle(points.begin(), points.end(), gen);

  cdt_.clear();

  QString fault_filename = filename;
  fault_filename.replace(QRegularExpression("\\.data$"), "_fault.data");
  std::ifstream ins2(qPrintable(fault_filename));
  if(!ins2) {
    return;
  }
  int num_lines = 0;
  ins2 >> num_lines;
  std::vector<int> num_vertex_per_line(num_lines);
  for(int n = 0; n < num_lines; ++n) {
    ins2 >> num_vertex_per_line[n];
  }

  CGAL::Bbox_2 b;
  for(int i = 0; i < num_lines; ++i) {
    Point_2 p, q;
    ins2 >> p;
    if(i == 0) {
      b = p.bbox();
    } else {
      b = b + p.bbox();
    }
    for(int j = 1; j < num_vertex_per_line[i]; ++j) {
      ins2 >> q;
      cdt_.insert_constraint(p, q);
      p = q;
      b = b + p.bbox();
    }
  }

  for(const auto& point : points) {
    if(CGAL::do_overlap(b, point.bbox())) {
      cdt_.insert(point);
    }
  }

  xmax = b.xmax();
  xmin = b.xmin();
  ymax = b.ymax();
  ymin = b.ymin();

  dx = (xmax - xmin) / 20.0;
  dy = (ymax - ymin) / 20.0;
  xmin -= dx;
  ymin -= dy;
  xmax += dx;
  ymax += dy;
  Point_2 bl(xmin, ymin);
  Point_2 br(xmax, ymin);
  Point_2 tl(xmin, ymax);
  Point_2 tr(xmax, ymax);
  cdt_.insert_constraint(bl, br);
  cdt_.insert_constraint(br, tr);
  cdt_.insert_constraint(tr, tl);
  cdt_.insert_constraint(tl, bl);
}

void MainWindow::updateStepLength(int value) {
  step_length_ = value;
}

void MainWindow::updateTimerInterval(int value) {
  timer_interval_ = value;
  if(timer_->isActive()) {
    timer_->start(timer_interval_);
  }
}

void MainWindow::updateCriteria() {
  criteria_.set_bound(angle_bound_box_->value());
  criteria_.set_size_bound(size_bound_box_->value());
  criteria_.set_local_size(under_mouse_box_->isChecked());
  if(mesher_ != nullptr) {
    mesher_->set_criteria(criteria_);
  }
}

void MainWindow::onUnderMouseToggled(bool checked) {
  if(checked) {
    if(mesher_ == nullptr) {
      mesher_ = createMesher();
      mesher_->init();
      updateMesherState();
    }
    criteria_.set_local_size(true);
    mesher_->set_criteria(criteria_);
    scene_.installEventFilter(local_refiner_);
  } else {
    criteria_.set_local_size(false);
    if(mesher_ != nullptr) {
      mesher_->set_criteria(criteria_);
    }
    scene_.removeEventFilter(local_refiner_);
    local_refiner_->reset();
  }
}

void MainWindow::applyDisplaySettings() {
  dgi_->setEdgesPen(pen(display_settings_->edgeColor(), display_settings_->edgeWidth()));
  dgi_->setConstraintsPen(pen(display_settings_->constraintColor(), display_settings_->constraintWidth()));
  dgi_->setVerticesPen(pen(display_settings_->vertexColor(), display_settings_->vertexWidth()));
  dgi_->setFacesInDomainBrush(QBrush(display_settings_->facesColor()));
  dgi_->setSeedsPen(pen(display_settings_->seedsColor(), display_settings_->seedsWidth()));
  bad_faces_item_->setBrush(QBrush(display_settings_->badFacesColor()));
  encroached_edges_item_->setPen(
      pen(display_settings_->encroachedEdgesColor(), display_settings_->encroachedEdgesWidth()));
  clusters_item_->setPen(pen(display_settings_->clustersColor(), display_settings_->clustersWidth()));
  clusters_item_->setBrush(QBrush(display_settings_->clustersColor()));

  action_show_edges_->setChecked(display_settings_->showEdges());
  action_show_constraints_->setChecked(display_settings_->showConstraints());
  action_show_vertices_->setChecked(display_settings_->showVertices());
  action_show_faces_->setChecked(display_settings_->showFaces());
  action_show_seeds_->setChecked(display_settings_->showSeeds());
  action_show_bad_faces_->setChecked(display_settings_->showBadFaces());
  action_show_encroached_edges_->setChecked(display_settings_->showEncroachedEdges());
  action_show_clusters_->setChecked(display_settings_->showClusters());

  display_settings_->saveSettings();
  dgi_->modelChanged();
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Mesh_2 demo");

  CGAL_QT_INIT_RESOURCES;

  MainWindow window;
  window.show();

  for(const auto& arg : app.arguments().sliced(1)) {
    window.openTriangulation(arg);
  }

  return app.exec();
}

#include "mesh_2_demo.moc"
