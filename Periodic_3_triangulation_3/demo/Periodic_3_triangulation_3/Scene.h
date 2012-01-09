#ifndef SCENE_H
#define SCENE_H

#include "ui_MainWindow.h"
#include "Scene_utils.h"
#include <fstream>
#include <QObject>
#include <QFileDialog>
#include <CGAL/glu.h>

class Scene : public QObject
{
  Q_OBJECT

private:
  typedef qglviewer::Vec Vec;
  typedef std::set<Segment, Compare_segment<Segment> > Segment_set;

  enum Init {EMPTY, GRID, SINGLE, PLANE, RANDOM};

  enum Draw_type { SEGMENT, TRIANGLE, TETRAHEDRON };

  enum Primitive_colors {
    VERTEX_COLOR,
    EDGE_COLOR,
    DOMAIN_COLOR,
    FLYING_BALL_COLOR,
    CLIPPING_COLOR
  };

public:
  Scene(Ui::MainWindow* ui_) : ui(ui_), p3dt(),
				  moving_point(Point(0.2,0.2,0.2)) {

    flying_ball = ui->actionFlying_ball->isChecked();

    dlocate = ui->actionPoint_location->isChecked();
    dconflict = ui->actionConflict_region->isChecked();

    wireframe = ui->actionWireframe->isChecked();
    in_plane = ui->actionPlanar_triangulation->isChecked();
    bool c1 = ui->actionDraw_1_sheeted_covering->isChecked();
    bool cd = ui->actionDraw_bordering_cells_multiply->isChecked();

    if ( c1 &&  cd) it_type = P3DT::UNIQUE_COVER_DOMAIN;
    if ( c1 && !cd) it_type = P3DT::UNIQUE;
    if (!c1 &&  cd) it_type = P3DT::STORED_COVER_DOMAIN;
    if (!c1 && !cd) it_type = P3DT::STORED;

    if ( ui->actionDraw_segments->isChecked() )   draw_type = SEGMENT;
    if ( ui->actionDraw_triangles->isChecked() )  draw_type = TRIANGLE;
    if ( ui->actionDraw_tetrahedra->isChecked() ) draw_type = TETRAHEDRON;

    ddomain = ui->actionDraw_cube_square->isChecked();
    cube_clipping = ui->actionClip_along_the_cube_square->isChecked();
    two_color_clipping = ui->action2_color_clipping->isChecked();

    materials[VERTEX_COLOR]      = "Gold";
    materials[EDGE_COLOR]        = "Green";
    materials[DOMAIN_COLOR]      = "Black plastic";
    materials[FLYING_BALL_COLOR] = "Red";
    materials[CLIPPING_COLOR]    = "Silver";

    // Timer for flying ball
    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(update_position()));
    timer->start(100);
  }

  ~Scene() {
    gluDeleteQuadric(pQuadric);
  }

public slots:
  void init();
  void draw();

  void load_points(const QString& fileName);
  void update_position();

  void init_scene_empty() { init_scene(EMPTY); }
  void init_scene_single() { init_scene(SINGLE); }
  void init_scene_random() { init_scene(RANDOM); }
  void init_scene_plane() { init_scene(PLANE); }
  void init_scene_grid() { init_scene(GRID); }
  void toggle_flying_ball(bool on) {
    flying_ball = on;
    ui->viewer->update();
  }
  void toggle_pause(bool on) {
    if (on) timer->stop();
    else timer->start();
  }

  void insert_mp() {
    insert_point(moving_point);
    QString str;
    ui->viewer->displayMessage(str.sprintf("Added point (%f, %f, %f)",
		   moving_point.x(),moving_point.y(),moving_point.z()));
  }

  void insert_random() {
    RandPts rp(0.5);
    Point pt = *rp+Vector(0.5,0.5,0.5);
    rp++;
    insert_point(Point(pt.x(),pt.y(),(in_plane? 0.0:pt.z())));
    QString str;
    ui->viewer->displayMessage(str.sprintf("Added point (%f, %f, %f)",
		   pt.x(),pt.y(),(in_plane? 0.0:pt.z())));
  }
  void insert_point(Point p) {
    bool temp_flags[] = {dlocate, dconflict};
    dlocate = dconflict = false;
    make_draw_list();
    p3dt.insert(p);
    dlocate = temp_flags[0];
    dconflict = temp_flags[1];
    make_draw_list();
  }
  void grab_image() {
    ui->viewer->openSnapshotFormatDialog();
    ui->viewer->saveSnapshot(false);
  }
  void toggle_dlocate(bool on) { 
    dlocate = on;
    ui->viewer->update();
  }
  void toggle_dconflict(bool on) {
    dconflict = on; 
    ui->viewer->update();
  }
  void toggle_wireframe(bool on) {
    wireframe = on;
    ( on ? glDisable(GL_LIGHTING) : glEnable(GL_LIGHTING) );
    gl_draw_domain();
    make_draw_list();
  }
  void toggle_in_plane(bool on) {
    in_plane = on;
    gl_draw_domain();
    make_draw_list();
  }
  // TODO: find radio button functionality within menus in QtDesigner
  void toggle_force_1cover(bool on) {
    if (ui->actionDraw_bordering_cells_multiply->isChecked())
      it_type = ( on ? P3DT::UNIQUE_COVER_DOMAIN : P3DT::STORED_COVER_DOMAIN );
    else
      it_type = ( on ? P3DT::UNIQUE : P3DT::STORED);
    make_draw_list();
  }
  void toggle_multiple_cells(bool on) {
    if (ui->actionDraw_1_sheeted_covering->isChecked())
      it_type = ( on ? P3DT::UNIQUE_COVER_DOMAIN : P3DT::UNIQUE );
    else
      it_type = ( on ? P3DT::STORED_COVER_DOMAIN : P3DT::STORED );
    make_draw_list();
  }
  void trigger_draw_type_segment() {
    ui->actionDraw_segments->setChecked(true);
    ui->actionDraw_triangles->setChecked(false);
    ui->actionDraw_tetrahedra->setChecked(false);
    draw_type = SEGMENT;
    make_draw_list();
  }
  void trigger_draw_type_triangle() {
    ui->actionDraw_segments->setChecked(false);
    ui->actionDraw_triangles->setChecked(true);
    ui->actionDraw_tetrahedra->setChecked(false);
    draw_type = TRIANGLE;
    make_draw_list();
  }
  void trigger_draw_type_tetrahedron() {
    ui->actionDraw_segments->setChecked(false);
    ui->actionDraw_triangles->setChecked(false);
    ui->actionDraw_tetrahedra->setChecked(true);
    draw_type = TETRAHEDRON;
    make_draw_list();
  }
  void toggle_ddomain(bool on) { ddomain = on; }
  void toggle_cube_clipping(bool on) {
    ui->action2_color_clipping->setEnabled(on);
    if (!on) ui->action2_color_clipping->setChecked(false);
    cube_clipping = on;
    two_color_clipping = false;
    if (on && ui->action2_color_clipping->isChecked())
      toggle_two_color_clipping(true);
    else
      make_draw_list();
  }
  void toggle_two_color_clipping(bool on) {
    two_color_clipping = on;
    make_draw_list();
  }
  bool load_points() { 
    QString fileName = QFileDialog
      ::getOpenFileName(ui->centralWidget, tr("Open point set"),
			".", tr("All files (*)"));
    if(! fileName.isEmpty()){
      load_points(fileName);
    }
   return true;
  }

signals:
  void message(const QString & message, int timeout = 0 ); 
  void loaded_points(const QString &n);

private:
  // Scene management helpers
  void make_draw_list();
  void init_scene(Init sceneID);
  void gl_draw_square();
  void gl_draw_cube();
  void gl_draw_domain() {
    if (in_plane) gl_draw_square();
    else gl_draw_cube();
  }

  // Helper functions
  void get_tri_offsets(const Cell_handle ch, int i,
		       Offset& off0, Offset& off1, Offset& off2) const;
  void get_tet_offsets(const Cell_handle ch, Offset& off0, Offset& off1,
		       Offset& off2, Offset& off3) const;
  int get_tri_drawing_offsets(const Cell_handle ch, int i) const;
  int get_tet_drawing_offsets(const Cell_handle ch) const;
  Tetrahedron construct_tetrahedron(const Cell_handle ch, const Offset& off0,
      const Offset& off1,const Offset& off2,const Offset& off3, int off) const;
  Triangle construct_triangle(const Cell_handle ch, int i, 
      const Offset& off0,const Offset& off1,const Offset& off2, int off) const;

  // Drawing functions
  void gl_draw_vertex(const Point& p, const FT& radius) const;
  void gl_draw_edge(Point p1, Point p2, FT radius );

  void primitives_from_geom_it(Segment_set& sset);
  void segment_clipping(Segment_set& sset);
  void segment_2color_clipping (Segment_set& sset);

  void gl_draw_location();
  void gl_draw_conflict();

  void change_material(const QString& string);

private:
  Ui::MainWindow * ui;
  P3DT p3dt;
  Point moving_point;

  QString materials[6];
  QTimer* timer;
  GLuint l_triangulation, l_domain;
  GLUquadricObj* pQuadric;

  bool flying_ball;
  bool dlocate, dconflict;
  bool wireframe, in_plane;
  Iterator_type it_type;
  Draw_type draw_type;
  bool ddomain, cube_clipping, two_color_clipping;
};

#endif // SCENE_H
