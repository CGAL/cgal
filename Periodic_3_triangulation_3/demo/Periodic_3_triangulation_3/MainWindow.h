#include <QtGui/QMainWindow>
#include "Scene.h"

class MainWindow : public QMainWindow
{

public:
  MainWindow() {
    ui = new Ui::MainWindow;
    ui->setupUi(this);
    s = new Scene(ui);

    // QGLViewer drawing signals
    connect(ui->viewer, SIGNAL(viewerInitialized()), s, SLOT(init()));
    connect(ui->viewer, SIGNAL(drawNeeded()), s, SLOT(draw()));

    // divers
    connect(s, SIGNAL(message(const QString&, int)),
	    ui->statusBar, SLOT(showMessage(const QString&, int)));

    // File menu:
    connect(ui->actionLoad_Points, SIGNAL(triggered()),
	    s, SLOT(load_points()));
    connect(ui->actionExport_pov, SIGNAL(triggered()),
	    s, SLOT(export_pov()));

    // Init menu:
    connect(ui->actionEmpty_scene, SIGNAL(triggered()),
	    s, SLOT(init_scene_empty()));
    connect(ui->actionSingle_Point, SIGNAL(triggered()),
	    s, SLOT(init_scene_single()));
    connect(ui->actionRandom_Point_Set, SIGNAL(triggered()),
	    s, SLOT(init_scene_random()));
    connect(ui->actionRandom_Points_in_Plane, SIGNAL(triggered()),
	    s, SLOT(init_scene_plane()));
    connect(ui->actionPoint_grid, SIGNAL(triggered()),
	    s, SLOT(init_scene_grid()));

    // Actions menu:
    connect(ui->actionFlying_ball, SIGNAL(toggled(bool)),
	    s, SLOT(toggle_flying_ball(bool)));
    connect(ui->actionPause, SIGNAL(toggled(bool)),
	    s, SLOT(toggle_pause(bool)));

    connect(ui->actionInsert_point, SIGNAL(triggered()),
	    s, SLOT(insert_mp()));
    connect(ui->actionInsert_random_point, SIGNAL(triggered()),
	    s, SLOT(insert_random()));

    connect(ui->actionGrab_image, SIGNAL(triggered()),
	    s, SLOT(grab_image()));

    // Features menu:
    connect(ui->actionPoint_location, SIGNAL(toggled(bool)),
	    s, SLOT(toggle_dlocate(bool)));
    connect(ui->actionConflict_region, SIGNAL(toggled(bool)),
	    s, SLOT(toggle_dconflict(bool)));
    connect(ui->actionHole, SIGNAL(toggled(bool)),
	    s, SLOT(toggle_dhole(bool)));
    connect(ui->actionStar, SIGNAL(toggled(bool)),
	    s, SLOT(toggle_dstar(bool)));

    // Options menu:
    connect(ui->actionWireframe, SIGNAL(toggled(bool)),
	    s, SLOT(toggle_wireframe(bool)));
    connect(ui->actionPlanar_triangulation, SIGNAL(toggled(bool)),
	    s, SLOT(toggle_in_plane(bool)));

    connect(ui->actionDraw_1_sheeted_covering, SIGNAL(toggled(bool)),
	    s, SLOT(toggle_force_1cover(bool)));
    connect(ui->actionDraw_bordering_cells_multiply, SIGNAL(toggled(bool)),
	    s, SLOT(toggle_multiple_cells(bool)));

    connect(ui->actionDraw_segments, SIGNAL(triggered()),
	    s, SLOT(trigger_draw_type_segment()));
    connect(ui->actionDraw_triangles, SIGNAL(triggered()),
	    s, SLOT(trigger_draw_type_triangle()));
    connect(ui->actionDraw_tetrahedra, SIGNAL(triggered()),
	    s, SLOT(trigger_draw_type_tetrahedron()));

    connect(ui->actionDraw_cube_square, SIGNAL(toggled(bool)),
	    s, SLOT(toggle_ddomain(bool)));

    connect(ui->actionClip_along_the_cube_square, SIGNAL(toggled(bool)),
	    s, SLOT(toggle_cube_clipping(bool)));
    connect(ui->action2_color_clipping, SIGNAL(toggled(bool)),
	    s, SLOT(toggle_two_color_clipping(bool)));

  }

  ~MainWindow() {
    delete(ui);
    delete(s);
  }

  Ui::MainWindow* ui;
  Scene* s;
};

