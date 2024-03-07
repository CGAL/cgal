#include "window.h"

DemoWindowItem::DemoWindowItem() : CGAL::Qt::GraphicsItem(){
  // Clear
  _edges.clear();

  // Prepare the pens
  _poincare_disk_pen.setStyle(Qt::SolidLine);
  _poincare_disk_pen.setWidth(8);
  _poincare_disk_pen.setBrush(Qt::black);

  _edges_pen.setStyle(Qt::SolidLine);
  _edges_pen.setWidth(6);
  _edges_pen.setBrush(Qt::blue);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DemoWindowItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget){
  // 1. Draw the poincaré disk
  QRectF circle_rect = QRectF(-_poincare_disk_radius_in_pixels-3, -_poincare_disk_radius_in_pixels-3, 2*_poincare_disk_radius_in_pixels+6, 2*_poincare_disk_radius_in_pixels+6);
  painter->setPen(_poincare_disk_pen);
  painter->setBrush(QBrush());
  painter->drawEllipse(circle_rect);

  // 2. Draw the edges
  painter->setBrush(QBrush());
  painter->setPen(_edges_pen);
  for (int i=0; i<_edges.size(); i++){
    draw_edge(painter, _edges[i].first, _edges[i].second);
  }
}

QRectF DemoWindowItem::boundingRect() const{
  return QRectF(-_poincare_disk_radius_in_pixels-3, -_poincare_disk_radius_in_pixels-3, _poincare_disk_radius_in_pixels+6, _poincare_disk_radius_in_pixels+6);
}

void DemoWindowItem::modelChanged(){} // Only used by Qt : we don't need to fill it

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DemoWindowItem::draw_triangulation(Triangulation& triangulation){
  typedef std::vector<std::tuple<typename Triangulation::Combinatorial_map_with_cross_ratios::Dart_const_handle, typename Traits::Point_2,typename Traits::Point_2,typename Traits::Point_2>> RealizationVector;
  RealizationVector realized_triangles;
  realized_triangles = triangulation.lift();

  Complex point_1, point_2, point_3;
  for (typename RealizationVector::iterator it = realized_triangles.begin(); it != realized_triangles.end(); ++it){
    point_1 = std::get<1>(*it).get_z();
    point_2 = std::get<2>(*it).get_z();
    point_3 = std::get<3>(*it).get_z();

    _edges.push_back(std::make_pair(point_1, point_2));
    _edges.push_back(std::make_pair(point_2, point_3));
    _edges.push_back(std::make_pair(point_3, point_1));
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DemoWindowItem::draw_point(QPainter* painter, Complex position){
  // First convert the complex number of coordinates in doubles, well-scaled
  double point_x = _poincare_disk_radius_in_pixels * CGAL::to_double(position.real());
  double point_y = _poincare_disk_radius_in_pixels * CGAL::to_double(position.imag());

  // Then draw a small circle
  QRectF circle_rect = QRectF(point_x-1, point_y-1, 3, 3);
  painter->drawEllipse(circle_rect);

}

void DemoWindowItem::draw_edge(QPainter* painter, Complex source, Complex target){
  // First convert the points coordinates to doubles

  double src_x = CGAL::to_double(source.real());
  double src_y = CGAL::to_double(source.imag());

  double tar_x = CGAL::to_double(target.real());
  double tar_y = CGAL::to_double(target.imag());

  // 0. If src and tar are too colinear or too close from each other then draw a line
  double determinant = src_x*tar_y - src_y*tar_x; // determinant of the matrix whose columns are the vectors src and tar : indicates colinearity
  double distance_squared = (src_x-tar_x)*(src_x-tar_x) + (src_y-tar_y)*(src_y-tar_y);
  if ((std::abs(determinant) < computation_treshold_squared) || (distance_squared < computation_treshold_squared)){
      // src and tar are too colinear or too close from each other
      draw_line(painter, src_x, src_y, tar_x, tar_y);
      return;
  }

  // 1. Compute the center of the circle supporting the geodesic between src and tar

  // 1.a Inverse src and tar with respect to the unit circle and find the euclidean midpoints of the segments between respectively
  //       src and it's inversion, and tar and it's inversion

  double src_norm_2 = src_x*src_x + src_y*src_y; // Can't be too close to zero because determinant was not
  double tar_norm_2 = tar_x*tar_x + tar_y*tar_y; // Can't be too close to zero because determinant was not

  double src_inv_x = src_x / src_norm_2;
  double src_inv_y = src_y / src_norm_2;
  double tar_inv_x = tar_x / tar_norm_2;
  double tar_inv_y = tar_y / tar_norm_2;

  // coordinates of the euclidean midpoints of the segments [src, src_inv] and [tar, tar_inv]
  double src_mid_x = (src_x + src_inv_x) / 2;
  double src_mid_y = (src_y + src_inv_y) / 2;
  double tar_mid_x = (tar_x + tar_inv_x) / 2;
  double tar_mid_y = (tar_y + tar_inv_y) / 2;

  // 1.b Solve a system to find the intersection (center_x, center_y) of the bisectors of the two segments [src, src_inv] and [tar, tar_inv] :
  //      (center_x \\ center y) = (a & b \\ c & d)^{-1} \times (u_x \\ u_y)


  // 1.b.i define the system

  double a = src_x;
  double b = src_y;
  double c = tar_x;
  double d = tar_y;

  double u_x = a*src_mid_x + b*src_mid_y;
  double u_y = c*tar_mid_x + d*tar_mid_y;

  // 1.b.ii solve the system (just a matrix inversion)
  double det = a*d-b*c; // Can't be too close to zero...
  double center_x = (d*u_x - b*u_y) / det;
  double center_y = (-c*u_x + a*u_y) / det;

  // 2. draw the arc supported by the circle whose center is (center_x, center_y) and whose extremities are src and tar
  draw_arc(painter, src_x, src_y, tar_x, tar_y, center_x, center_y);
}


void DemoWindowItem::draw_line(QPainter* painter, double point_1_x, double point_1_y, double point_2_x, double point_2_y){

  // Convert to doubles and scale by the radius of the poincaré disk
  double src_x = _poincare_disk_radius_in_pixels * point_1_x;
  double src_y = _poincare_disk_radius_in_pixels * point_1_y;
  double tar_x = _poincare_disk_radius_in_pixels * point_2_x;
  double tar_y = _poincare_disk_radius_in_pixels * point_2_y;

  // Actual drawing
  QLineF line (src_x, src_y, tar_x, tar_y);
  painter->drawLine(line);
}


void DemoWindowItem::draw_arc(QPainter* painter, double point_1_x, double point_1_y,
        double point_2_x, double point_2_y, double center_x, double center_y){
  // Draws the arc supported by the circle whose center is (center_x, center_y) and whose extremities are src and tar

  // 1. Scale by the radius of the poincaré disk

  double src_x = _poincare_disk_radius_in_pixels * point_1_x;
  double src_y = _poincare_disk_radius_in_pixels * point_1_y;

  double tar_x = _poincare_disk_radius_in_pixels * point_2_x;
  double tar_y = _poincare_disk_radius_in_pixels * point_2_y;

  double xc = _poincare_disk_radius_in_pixels * center_x;
  double yc = _poincare_disk_radius_in_pixels * center_y;

  // 2. Define the radius of the circle and the box [xm, xM] \times [ym, yM] bounding the circle

  double circle_radius = sqrt((point_1_x-center_x)*(point_1_x-center_x) + (point_1_y-center_y)*(point_1_y-center_y));

  double xm = _poincare_disk_radius_in_pixels * (center_x - circle_radius);
  double xM = _poincare_disk_radius_in_pixels * (center_x + circle_radius);
  double ym = _poincare_disk_radius_in_pixels * (center_y - circle_radius);
  double yM = _poincare_disk_radius_in_pixels * (center_y + circle_radius);

  // If the source and the target are too close from each other (less than 10 pixels) or if the circle is very big then just draw a line
  double dist_sq = (src_x-tar_x)*(src_x-tar_x) + (src_y - tar_y)*(src_y-tar_y);
  double rad_sq = (xM-xc)*(xM-xc) + (yM-yc)*(yM-yc);
  if ( (dist_sq < 100) ||  (rad_sq > 1000 * dist_sq) ){
    QLineF line (src_x, src_y, tar_x, tar_y);
    painter->drawLine(line);
    return;
  }

  // 3. Compute angles (needed because we will draw using QPainter::drawArc)

  // src_angle is the argument, in degrees, of point_1 - center
  // tar_angle is the argument, in degrees, of point_2 - center
  double src_angle = deg_angle(src_x - xc, src_y - yc);
  double tar_angle = deg_angle(tar_x - xc, tar_y - yc);
  src_angle = 360 - src_angle; // Because of y-axis inversion
  tar_angle = 360 - tar_angle; // Because of y-axis-inversion

  // Compute the sweep angle (see QPainter::drawArc)
  double sweep_angle = tar_angle - src_angle;
  while (sweep_angle > 180)
    sweep_angle -= 360;
  while (sweep_angle < -180)
    sweep_angle += 360;

  // 4. Actual Drawing

  QRectF bbox_rect (xm, ym, xM-xm, yM-ym);
  painter->drawArc(bbox_rect, src_angle*16, sweep_angle*16);
}

double DemoWindowItem::deg_angle(double x, double y){
    // To avoid problems when further division by x (ok since x^2 + y^2 not too small) :
    if (x*x <  computation_treshold_squared){
      if (y>0) return 90;
      return -90;
    }

    double angle = 180. * std::atan(y / x) / M_PI;
    if (x < 0){
      return angle + 180.;
    }
    return angle;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DemoWindow::DemoWindow() : DemosMainWindow(){
  setupUi(this); // Method automatically generated by the ui file and here inherited from Ui::MainWindow. Builds the window and the contents for us...
  this->graphicsView->setScene(&_scene); // ... in particular graphicsView is already constructed : we just put our scene in it and then do things within the scene
  _scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  _scene.setSceneRect(-600, -600, 1200, 1200);
  _item = new DemoWindowItem();
  _scene.addItem(_item);
  this->graphicsView->scale(0.5, -0.5); // Y-axis inversion

  setWindowTitle("Hyperbolic surfaces triangulation 2 Demo");
}

DemoWindowItem& DemoWindow::item(){
  return *_item;
}

void DemoWindow::keyPressEvent(QKeyEvent* event) {}
