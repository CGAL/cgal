// ============================================================================
//
// Copyright (c) 2001-2010 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/).
//
// ----------------------------------------------------------------------------
//
// Library       : AlciX
// File          : demos/xalci/misc.cpp
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.12 $
// Revision_date : $Date: 2009-06-30 13:14:58 $
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-inf.mpg.de>
//
// ============================================================================

/*!\brief
 * declares miscellaneous routines in a separate file (to speed-up
 * compilation)
 */

#define NDEBUG 1

//#define CGAL_NO_LEDA
//#define  CGAL_POLYNOMIAL_USE_NTL_MUL
#define CGAL_ACK_CHECK_SQUAREFREENESS 0

#include <CGAL/basic.h>
#include <CGAL/Timer.h>

#include <set>
#include <ctime>
#include <fstream>

#define CGAL_CKvA_LINE_THICKNESS 4
#define CGAL_CKvA_NO_AXES 1

#include "include/misc.h"
#include "include/xalci.h"

#include "axis.xpm"

#include <CGAL/Curved_kernel_via_analysis_2/Qt_widget_Curve_renderer_2.h>

#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>

#include <CGAL/Polynomial_type_generator.h>

CGAL::Timer refine_timer;

QColor rasterize_colors[] = {
    QColor(139, 0, 0),
    QColor(138, 43, 226),
    QColor(95, 158, 160),
    QColor(0, 0, 139),
    QColor(205, 149, 12),
    QColor(0, 100, 0),
    QColor(139, 0, 139),
    QColor(85, 107, 47),
    QColor(255, 127, 0),
    QColor(0, 206, 209),
    QColor(238, 18, 137),
    QColor(238, 99, 99),
    QColor(205, 55, 0)
};
int n_rast_colors = 13;

static CGAL::Bbox_2 bbox(0.0, 0.0, 0.0, 0.0);
//static SoX::Subdiv_renderer subdiv_renderer;

Object_vector cad_objects, oc_objects, arr_objects;
Object_vector* curr_objects;

Graphic_layer *subdiv_layer;
Layers cad_layers, oc_layers, arr_layers;

QPixmap *subdiv_plot, *arcs_plot;
bool subdiv_layer_changed = true;

CGAL::Timer timer;

void Graphic_layer::draw()
{
    QPainter *ppnt = &widget->get_painter();
    QPen old_pen = ppnt->pen();

    if(index == -2) {
        CGAL::Bbox_2 new_box(widget->x_min(), widget->y_min(),
                             widget->x_max(), widget->y_max());
        if(bbox != new_box) {
            bbox = new_box;

            //subdiv_renderer.setup(bbox,
                // widget->width(), widget->height());
            subdiv_layer_changed = true;
        }

        if(subdiv_layer_changed) {
            QPainter offscreen(subdiv_plot);
            //offscreen.setPen(QPen(Qt::black, 2, Qt::SolidLine, Qt::SquareCap,
              //  Qt::MiterJoin));
            subdiv_plot->fill();

            subdiv_layer_changed = false;
            typedef std::pair<int, int> Int_pair;
            std::list<Int_pair> points;

            //subdiv_renderer.draw(std::back_inserter(points));
            if(!points.empty()) {

                offscreen.setPen(QPen(Qt::black ,1));
                std::list<Int_pair>::iterator it = points.begin();
                while(it != points.end()) {
                    offscreen.drawPoint(it->first, it->second);
                    it++;
                }
            }
            subdiv_layer_changed = false;
        }
        ppnt->drawPixmap(0,0,*subdiv_plot);

    } else if(index == -1) { // this layer is dedicated to axis drawing

#if !CGAL_CKvA_NO_AXES
        RasterOp old_raster = widget->rasterOp();
        widget->setRasterOp(XorROP);
        ppnt->setPen(QPen(QColor(150,150,0),1,Qt::DashDotDotLine));
        ppnt->moveTo(0,widget->y_pixel(0));
        ppnt->lineTo(widget->width(),widget->y_pixel(0));
        ppnt->moveTo(widget->x_pixel(0),0);
        ppnt->lineTo(widget->x_pixel(0),widget->height());
        widget->setRasterOp(old_raster);
#endif

    } else {
        widget->setRasterOp(CopyROP);
        ppnt->setPen(QPen(rasterize_colors[color_index % n_rast_colors],
             CGAL_CKvA_LINE_THICKNESS,  Qt::SolidLine, Qt::RoundCap,
                 Qt::MiterJoin));

        const CGAL::Object obj = (*curr_objects)[index];
        Arc_2 arc;
        Point_2 pt;

//#if !XALCI_USE_FLAT_COLOR_SCHEME
         QBrush b1(Qt::black);
/*#else
        QBrush b1(Qt::NoBrush);
#endif*/
        ppnt->setBrush(b1);

        timer.start();
        if(CGAL::assign(arc, obj))
           *widget << arc;
        else if(CGAL::assign(pt, obj))
           *widget << pt;
        else
            CGAL_error_msg("Malformed object found..\n");
        timer.stop();
    }
    ppnt->setPen(old_pen);
}

void xAlci_main_window::arr_activate_layers() {

    int i = 0;
// int n_nodes =  arr_node_list->numRows();

    for(Layers::iterator it = arr_layers.begin(); it != arr_layers.end();
            it++, i++) {
        /*if(i < n_nodes) {
            if(arr_node_list->isSelected(i)||arr_complete_check->isChecked())
                (*it)->activate();
            else if((*it)->is_active())
                (*it)->deactivate();
        } else {*/
        if(arr_edge_list->isSelected(i) || arr_complete_check->isChecked())
            (*it)->activate();
        else if((*it)->is_active())
            (*it)->deactivate();
    }
}


void xAlci_main_window::cad_activate_layers()
{
    Layers::iterator it = cad_layers.begin();
    int i = 0;
    while(it != cad_layers.end()) {
        if(cad_seg_list->isSelected(i)||cad_complete_check->isChecked()) {
            (*it)->activate();
        } else if((*it)->is_active()) {
            (*it)->deactivate();
        }
        it++; i++;
    }
}

void xAlci_main_window::oc_activate_layers()
{
    Layers::iterator it = oc_layers.begin();
    int i = 0;
    while(it != oc_layers.end()) {
        if(oc_seg_list->isSelected(i) || oc_complete_check->isChecked()) {
            (*it)->activate();
        } else if((*it)->is_active()) {
            (*it)->deactivate();
        }
        it++; i++;
    }
}

void xAlci_main_window::arr_rasterize_click() {
  arr_activate_layers();

  timer.reset();
  widget->redraw();
  std::cout << "\n\nRasterize elapsed time: " << timer.time() << std::endl;
}

void xAlci_main_window::cad_rasterize_click()
{
  cad_activate_layers();

  timer.reset();
  widget->redraw();
  std::cout << "\n\nRasterize elapsed time: " << timer.time() << std::endl;
}

void xAlci_main_window::oc_rasterize_click()
{
    if(cur_method == 0) {
        oc_activate_layers();
    } else {
//         CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
        ::CGAL::IO::set_pretty_mode(std::cout);
        Poly_int2 f;

        if(!input_poly(f, oc_input->text().ascii()))
            return;

        //subdiv_renderer.set_polynomial(f);
        subdiv_layer_changed = true;
    }

    refine_timer.reset();
    timer.reset();

    widget->redraw();
    refine_timer.stop();

    std::cout << "\n\nrefine end-points: " << refine_timer.time() <<
        std::endl;
    std::cout << "\n\nRasterize elapsed time: " << timer.time() << std::endl;
}

bool xAlci_main_window::input_poly(Poly_int2& p, const char *ascii) {

    if(ascii == NULL)
        return false;

    typedef CGAL::Polynomial_type_generator< Rational, 2 >::Type Poly_rat_2;

    CGAL::Polynomial_parser_d< Poly_rat_2, CGAL::Mixed_floating_point_parser_policy< Poly_rat_2 > > parser;
    std::string str(ascii);

    Poly_rat_2 prat;
    if(parser(str, prat)) {

        typedef CGAL::Fraction_traits< Poly_rat_2 > FTraits;
        FTraits::Denominator_type det(1);
        FTraits::Decompose decompose;
        decompose(prat, p, det);

    } else {
        std::cerr << "Parser error, trying another format..\n";
        try {
            std::stringstream ss(str);
            ss >> p;
        }
        catch(...) {
            std::cerr << "Invalid format of polynomial..\n";
            return false;
        }
    }
#if !AcX_SQRT_EXTENSION && CGAL_ACK_CHECK_SQUAREFREENESS
    if(!CGAL::CGALi::is_square_free(p)) {
        p = make_square_free(p);
        std::cout << "squarefree part: " << p << std::endl;
    }
#endif
    return true;
}

void xAlci_main_window::print_arc(const Arc_2& arc, std::ostream& os) {

    print_endpoint(arc, CGAL::ARR_MIN_END, os);
    os << " - ";
    print_endpoint(arc, CGAL::ARR_MAX_END, os);

    if(arc.is_vertical())
        os << "; vertical";
    else
        os << "; arcno: " << arc.arcno();
}

void xAlci_main_window::print_endpoint(const Arc_2& arc,
    CGAL::Arr_curve_end end, std::ostream& os) {

    CGAL::Arr_parameter_space loc = arc.location(end);

    os << "(x: ";
    if(loc == CGAL::ARR_LEFT_BOUNDARY)
        os << "-oo";
    else if(loc == CGAL::ARR_RIGHT_BOUNDARY)
        os << "+oo";
    else
        os << CGAL::to_double(arc.curve_end_x(end));

    switch(loc) {
    case CGAL::ARR_BOTTOM_BOUNDARY:
        os << "; y: -oo)";
        break;

    case CGAL::ARR_TOP_BOUNDARY:
        os << "; y: +oo)";
        break;

    case CGAL::ARR_LEFT_BOUNDARY:
    case CGAL::ARR_RIGHT_BOUNDARY: {

        CGAL::Object obj =
            arc.curve().asymptotic_value_of_arc(loc, arc.arcno());

        CGAL::Arr_parameter_space loc2;
        if(CGAL::assign(loc2, obj))
            os << (loc2 == CGAL::ARR_BOTTOM_BOUNDARY ? "; y: -oo)" :
                "; y: +oo)");
        else {
            Kernel_2::Coordinate_1 y;
            if(CGAL::assign(y, obj))
                os << "; y: " << CGAL::to_double(y) << " (asym))";
            else
                CGAL_error_msg("Ill-typed object returned..\n");
        }
        break;
    }
    default:
        os << "; arcno: " << arc.curve_end(end).arcno() << ')';
    }
}

void xAlci_main_window::print_point(const Point_2& pt,
    std::ostream& os) {

    os << "(x: " << CGAL::to_double(pt.x()) <<
        "; arcno: " << pt.arcno() << ")";
}

void xAlci_main_window::arr_deactivate_layers()
{
    Layers::iterator it = arr_layers.begin();
    while(it != arr_layers.end())
        (*it++)->deactivate();
}

void xAlci_main_window::cad_deactivate_layers()
{
    Layers::iterator it = cad_layers.begin();
    while(it != cad_layers.end())
        (*it++)->deactivate();
}

void xAlci_main_window::oc_deactivate_layers()
{
    Layers::iterator it = oc_layers.begin();
    while(it != oc_layers.end())
        (*it++)->deactivate();
}

void xAlci_main_window::arr_complete_toggle(bool on) {
    if(on && !arr_objects.empty())
        arr_rasterize_btn->setEnabled(true);
}

void xAlci_main_window::cad_complete_toggle(bool on) {
    if(on && !cad_objects.empty())
        cad_rasterize_btn->setEnabled(true);
}

void xAlci_main_window::oc_complete_toggle(bool on) {
    if(on && !oc_objects.empty())
        oc_rasterize_btn->setEnabled(true);
}


void xAlci_main_window::cad_file_search_click() {
    QFileDialog file_dialog("", QString::null, central_widget, 0, true);

    if(file_dialog.exec() == QDialog::Accepted)
        cad_input->setText(file_dialog.selectedFile());
}


void xAlci_main_window::arr_file_search_click() {
    QFileDialog file_dialog("", QString::null, central_widget, 0, true);

    if(file_dialog.exec() == QDialog::Accepted)
        arr_input->setText(file_dialog.selectedFile());
}

void xAlci_main_window::arr_node_list_click() {
  arr_rasterize_btn->setEnabled(true);
}

void xAlci_main_window::arr_edge_list_click() {
  arr_rasterize_btn->setEnabled(true);
}

void xAlci_main_window::cad_seg_list_click()
{
    cad_rasterize_btn->setEnabled(true);
}

void xAlci_main_window::oc_seg_list_click()
{
    oc_rasterize_btn->setEnabled(true);
}

void xAlci_main_window::axis_toggle()
{
std::cerr << "im here\n";
    widget->lock();
    axis->draw();
    widget->unlock();
    widget->repaint(false);
}

void xAlci_main_window::setup(int w, int h)
{
    err_msg_dialog = new QErrorMessage(this);

    central_widget = new QWidget(this);
    setCentralWidget(central_widget);

    subdiv_plot = new QPixmap(w, h);
    arcs_plot = new QPixmap(w, h);
    subdiv_plot->fill();
    arcs_plot->fill();

    subdiv_layer = NULL;
    QBoxLayout *hbox = new QHBoxLayout(central_widget, 10, 10);
    widget = new CGAL::Qt_widget(central_widget);
    hbox->addWidget(widget,8);
        //QBoxLayout *vbox = new QVBoxLayout(0,0,5);
        //hbox->addLayout(vbox);

    tab_widget = new QTabWidget(central_widget);
    hbox->addWidget(tab_widget,4);
    one_curve_tab = new QFrame(tab_widget,"one_curve");
    cad_tab = new QFrame(tab_widget,"cad");
    arr_tab = new QFrame(tab_widget,"arrangement");

    curr_objects=&oc_objects;
    tab_widget->addTab(one_curve_tab,"one_curve");
    tab_widget->addTab(cad_tab,"cad");
    tab_widget->addTab(arr_tab,"arrangement");

    *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor(CGAL::white());
    resize(w,h);
    double ratio = 1.0;//(double)h/w;
    widget->set_window(-1, 1, -ratio, ratio, true);
        widget->setMouseTracking(TRUE);
    subdiv_layer = new Graphic_layer(widget, -2, 0);
    subdiv_layer->deactivate();
    axis = new Graphic_layer(widget, -1, 0);

        // LAYOUT CAD-TAB

    QBoxLayout* cad_vbox = new QVBoxLayout(cad_tab,10,10);
    cad_vbox->addWidget(new QLabel("<b>Input file:</b>",cad_tab));
    QHBox* cad_hbox1 = new QHBox(cad_tab);
    cad_input = new QLineEdit("", QString::null,cad_hbox1);

    cad_file_search = new QPushButton("Browse",cad_hbox1);
    cad_vbox->addWidget(cad_hbox1);
    QHBox* cad_hbox2 = new QHBox(cad_tab);

    cad_analyse_btn = new QPushButton("Analyse",cad_hbox2);
    cad_partial_selection = new QPushButton("Choose polynomials",cad_hbox2);

    cad_vbox->addWidget(cad_hbox2);
    cad_vbox->addWidget(new QLabel("<b>Curve segments:</b>",cad_tab));

    cad_seg_list = new QListBox(cad_tab);
    cad_seg_list->setSelectionMode(QListBox::Multi);
    cad_vbox->addWidget(cad_seg_list,6);
    cad_vbox->addWidget(new QLabel("<b>Curves:</b>",cad_tab));

    cad_curve_list = new QListBox(cad_tab);
    cad_curve_list->setSelectionMode(QListBox::Multi);

    cad_vbox->addWidget(cad_curve_list,6);
    cad_complete_check = new QCheckBox("rasterize complete cad",cad_tab);
    cad_vbox->addWidget(cad_complete_check);

    cad_rasterize_btn = new QPushButton("Rasterize",cad_tab);
    cad_vbox->addWidget(cad_rasterize_btn);

        // ARRANGEMENT TAB
    QBoxLayout* arr_vbox = new QVBoxLayout(arr_tab,10,10);
    arr_vbox->addWidget(new QLabel("<b>Input file:</b>",arr_tab));
    QHBox* arr_hbox1 = new QHBox(arr_tab);

    arr_input = new QLineEdit("", QString::null,arr_hbox1);
    arr_file_search = new QPushButton("Browse",arr_hbox1);
    arr_vbox->addWidget(arr_hbox1);

    arr_method = new QHButtonGroup(arr_tab);
    arr_method->setTitle("Arrangement data structure");
    arr_leda = new QRadioButton(arr_method);
    arr_leda->setText("LEDA");
    arr_cgal = new QRadioButton(arr_method);
    arr_cgal->setText("CGAL");
    arr_method->setExclusive(true);

    arr_vbox->addWidget(arr_method);
    QHBox* arr_hbox2 = new QHBox(arr_tab);

    arr_analyse_btn = new QPushButton("Analyse",arr_hbox2);
    arr_partial_selection = new QPushButton("Choose polynomials",arr_hbox2);

    arr_vbox->addWidget(arr_hbox2);
    arr_node_label = new QLabel("<b>Nodes:</b>",arr_tab);
    arr_vbox->addWidget(arr_node_label);

    arr_node_list = new QListBox(arr_tab);
    arr_node_list->setSelectionMode(QListBox::Multi);

    arr_vbox->addWidget(arr_node_list,6);
    arr_edge_label = new QLabel("<b>Edges:</b>",arr_tab);
    arr_vbox->addWidget(arr_edge_label);

    arr_edge_list = new QListBox(arr_tab);
    arr_edge_list->setSelectionMode(QListBox::Multi);

    arr_vbox->addWidget(arr_edge_list,6);
    arr_complete_check =
        new QCheckBox("rasterize complete arrangement", arr_tab);
    //        complete_check->setChecked(true);
    arr_vbox->addWidget(arr_complete_check);

    arr_rasterize_btn = new QPushButton("Rasterize",arr_tab);
    arr_vbox->addWidget(arr_rasterize_btn);

        // ONE CURVE TAB
    QBoxLayout* oc_vbox = new QVBoxLayout(one_curve_tab,10,10);
    oc_vbox->addWidget(new QLabel("<b>Input polynomial:</b>",one_curve_tab));

    oc_input = new QTextEdit("", QString::null,one_curve_tab);
    oc_vbox->addWidget(oc_input,6);
    oc_analyse_btn = new QPushButton("Analyse",one_curve_tab);
    oc_vbox->addWidget(oc_analyse_btn);
    oc_vbox->addWidget(new QLabel("<b>Curve segments:</b>",one_curve_tab));

    oc_seg_list = new QListBox(one_curve_tab);
    oc_seg_list->setSelectionMode(QListBox::Multi);
    oc_vbox->addWidget(oc_seg_list,6);

    oc_complete_check =
        new QCheckBox("rasterize complete curve",one_curve_tab);
    oc_vbox->addWidget(oc_complete_check);

    oc_method_box = new QComboBox("Rasterization method", one_curve_tab);
    oc_method_box->insertItem("Segment Renderer");
    oc_method_box->insertItem("Space Subdivision");
    oc_method_box->setEditable(false);

    oc_vbox->addWidget(oc_method_box);
    oc_rasterize_btn = new QPushButton("Rasterize",one_curve_tab);
    oc_vbox->addWidget(oc_rasterize_btn);


    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ),
        CTRL+Key_Q );
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), 0);
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    QToolBar *layers_toolbar =
    new QToolBar("Tools", this, QMainWindow::DockTop, TRUE, "Tools");
    QToolButton *axis_button = new QToolButton(QPixmap(axis_xpm),
               "Show axis", 0, this, SLOT(axis_toggle()), layers_toolbar);
    axis_button->setToggleButton(true);
    axis_button->toggle();

    connect(tab_widget,SIGNAL(currentChanged(QWidget*)),
        SLOT(tab_changed(QWidget*)));

    //connect(widget, SIGNAL(rangesChanged()), SLOT(rasterize()));

    connect(cad_analyse_btn, SIGNAL(clicked()), SLOT(cad_analyse_click()));
    connect(cad_partial_selection,SIGNAL(clicked()),
         SLOT(cad_partial_selection_click()));

    connect(cad_rasterize_btn, SIGNAL(clicked()), SLOT(cad_rasterize_click()));
    connect(cad_file_search,SIGNAL(clicked()), SLOT(cad_file_search_click()));
    connect(cad_complete_check, SIGNAL(toggled(bool)),
        SLOT(cad_complete_toggle(bool)));

    connect(cad_seg_list, SIGNAL(selectionChanged()),
        SLOT(cad_seg_list_click()));
    connect(cad_curve_list, SIGNAL(selectionChanged()),
        SLOT(cad_curve_list_click()));

    connect(arr_analyse_btn, SIGNAL(clicked()), SLOT(arr_analyse_click()));
    connect(arr_partial_selection,SIGNAL(clicked()),
        SLOT(arr_partial_selection_click()));

    connect(arr_rasterize_btn, SIGNAL(clicked()), SLOT(arr_rasterize_click()));
    connect(arr_file_search,SIGNAL(clicked()), SLOT(arr_file_search_click()));

    connect(arr_complete_check, SIGNAL(toggled(bool)),
        SLOT(arr_complete_toggle(bool)));

    connect(arr_node_list, SIGNAL(selectionChanged()),
        SLOT(arr_node_list_click()));
    connect(arr_edge_list, SIGNAL(selectionChanged()),
        SLOT(arr_edge_list_click()));

    connect(oc_analyse_btn, SIGNAL(clicked()), SLOT(oc_analyse_click()));
    connect(oc_rasterize_btn, SIGNAL(clicked()), SLOT(oc_rasterize_click()));
    connect(oc_complete_check, SIGNAL(toggled(bool)),
        SLOT(oc_complete_toggle(bool)));
    connect(oc_seg_list, SIGNAL(selectionChanged()),
        SLOT(oc_seg_list_click()));
    connect(oc_method_box, SIGNAL(activated(int)), this,
        SLOT(oc_switch_method(int)));

    stoolbar = new CGAL::Qt_widget_standard_toolbar(widget, this, "ST");
    cad_rasterize_btn->setEnabled(false);
        oc_rasterize_btn->setEnabled(false);
    arr_rasterize_btn->setEnabled(false);

    CORE::CORE_init(2);
    CORE::setDefaultPrecision(70, CORE::extLong::getPosInfty());
    tab_widget->showPage(arr_tab);
}


#include "misc.moc"
