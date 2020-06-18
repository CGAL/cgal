// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/).
//
// ----------------------------------------------------------------------------
//
// Library       : AlciX
// File          : demos/xalci/include/xalci.h
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.37 $
// Revision_date : $Date: 2009-06-30 13:14:59 $
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef XALCI_H
#define XALCI_H

#include <iostream>
#include <vector>

#include <qplatinumstyle.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qerrormessage.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qtimer.h>
#include <qlistbox.h>
#include <qlayout.h>
#include <qpushbutton.h>
#include <qlabel.h>
#include <qtextedit.h>
#include <qcheckbox.h>
#include <qcombobox.h>
#include <qpixmap.h>
#include <qpainter.h>
#include <qtabwidget.h>
#include <qhbuttongroup.h>
#include <qradiobutton.h>
#include <qvbox.h>
#include <qfiledialog.h>
#include <qprogressdialog.h>

#ifndef CGAL_ACK_DEBUG_FLAG
#define CGAL_ACK_DEBUG_FLAG 0
#endif

#ifndef CGAL_ACK_DEBUG_PRINT
#define CGAL_ACK_DEBUG_PRINT std::cout
#endif

#ifndef CGAL_ACK_MAKE_SQUARE_FREE
#define CGAL_ACK_MAKE_SQUARE_FREE 0
#endif

// #define NDEBUG 1

// #define CGAL_ACK_BENCHMARK_RES
// #define CGAL_ACK_BITSTREAM_USES_E08_TREE 0

#include "misc.h"

#include <CGAL/Arithmetic_kernel.h>
// #include <CGAL/Polynomial/sturm_habicht_sequence.h>

#include <CGAL/Algebraic_kernel_d/flags.h>
#include <CGAL/Algebraic_curve_kernel_2_generator.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h>

#include <CGAL/Bbox_2.h>

#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>

typedef CGAL::CORE_arithmetic_kernel AK;
typedef AK::Rational Rational;
typedef AK::Integer Integer;

#if !AcX_SQRT_EXTENSION

typedef Integer Coefficient;

#else

typedef CGAL::Sqrt_extension<CGAL::Sqrt_extension<CGAL::Sqrt_extension<Integer,Integer>,CGAL::Sqrt_extension<Integer,Integer> >,CGAL::Sqrt_extension<CGAL::Sqrt_extension<Integer,Integer>,CGAL::Sqrt_extension<Integer,Integer> > > Coefficient;

//typedef CGAL::Sqrt_extension<Integer, Integer> Coefficient;
#endif

typedef CGAL::Algebraic_curve_kernel_2_generator<Coefficient>
      ::Algebraic_curve_kernel_with_qir_and_bitstream_2  Kernel_2;

// typedef Kernel_2::Coordinate_1 Coordinate_1;
// typedef Kernel_2::Coordinate_2 Coordinate_2;
typedef Kernel_2::Curve_analysis_2 Curve_analysis_2;

// types of supporting polynomials
typedef Curve_analysis_2::Polynomial_2 Poly_int2;
typedef Poly_int2::NT Poly_int1;

typedef CGAL::Curved_kernel_via_analysis_2<Kernel_2> CKvA_2;
typedef CKvA_2::X_monotone_curve_2 Arc_2;
typedef CKvA_2::Point_2 Point_2;

typedef std::vector<CGAL::Object> Object_vector;

//typedef SoX::Subdivision_1<NiX::Interval, AC2> Subdiv_renderer;

class xAlci_main_window : public QMainWindow
{
    Q_OBJECT

public:

    xAlci_main_window(int w, int h) : bbox(0.0, 0.0, 0.0, 0.0), cur_method(0) {
        setup(w, h);
    }

    void visualize();

public slots:

    void cad_analyse_click();
    void cad_rasterize_click();
    void cad_seg_list_click();
    void cad_curve_list_click();
    void cad_file_search_click();
    void cad_partial_selection_click();
    void cad_complete_toggle(bool);

    void arr_analyse_click();
    void arr_rasterize_click();
    void arr_edge_list_click();
    void arr_node_list_click();
    void arr_file_search_click();
    void arr_partial_selection_click();
    void arr_complete_toggle(bool);

    void oc_analyse_click();
    void oc_rasterize_click();
    void oc_seg_list_click();
    void oc_complete_toggle(bool);
    void oc_switch_method(int);

    void axis_toggle();

    void new_instance()
    {
        widget->lock();
        widget->clear_history();
        widget->set_window(-1.1, 1.1, -1.1, 1.1);
        // set the Visible Area to the Interval
        widget->unlock();
    }

protected slots:

    void about()
    {
        QMessageBox::about(this, "About",
            "This program demonstrates the use of Segment renderern"
            "and space subdivision to visualize algebraic curvesn");
    }
    void aboutQt()
    {
        QMessageBox::aboutQt(this, "About Qt");
    }
    void howto()
    {
        CGAL::Qt_help_window *help =
          new CGAL::Qt_help_window("help/index.html", ".", 0, "help viewer");
        help->resize(400, 400);
        help->setCaption("Demo HowTo");
        help->show();
    }
    void new_window()
    {
        xAlci_main_window *ed = new xAlci_main_window(500, 500);
        ed->setCaption("Layer");
        ed->widget->clear_history();
        ed->widget->set_window(-1.1, 1.1, -1.1, 1.1);
        ed->show();
    }

    void tab_changed(QWidget*);
    void setup(int, int);

    void print_arc(const Arc_2& arc, std::ostream& os);
    void print_endpoint(const Arc_2& arc,
        CGAL::Arr_curve_end end, std::ostream& os);
    void print_point(const Point_2& pt, std::ostream& os);

    void arr_activate_layers();
    void cad_activate_layers();
    void oc_activate_layers();

    void arr_deactivate_layers();
    void cad_deactivate_layers();
    void oc_deactivate_layers();


protected:

    Poly_int2 make_square_free(const Poly_int2& poly);

    bool input_poly(Poly_int2& p, const char *ascii);

    template<typename OutputIterator>
    bool read_polys_from_file(QString filename,OutputIterator out);

    void cad_to_segments();
    void arr_compute_arrangement();

    CGAL::Bbox_2 bbox;

    std::vector< Poly_int2 > cad_curves, arr_curves;

    Curve_selection_dialog* curve_selection_dialog;

    QErrorMessage* err_msg_dialog;

    QPushButton *cad_analyse_btn, *cad_rasterize_btn, *cad_file_search,
          *cad_partial_selection,
          *arr_analyse_btn, *arr_rasterize_btn, *arr_file_search,
          *arr_partial_selection,
          *oc_analyse_btn, *oc_rasterize_btn;

    QHButtonGroup *arr_method;
    QRadioButton *arr_cgal, *arr_leda;
    QCheckBox *cad_complete_check,*oc_complete_check, *arr_complete_check;
    QListBox *cad_seg_list, *cad_curve_list,*oc_seg_list, *cad_ps_curve_list,
          *arr_edge_list, *arr_node_list;

    QTabWidget *tab_widget;
    QFrame* one_curve_tab, *cad_tab, *arr_tab;
    QLabel *arr_node_label,*arr_edge_label;
    std::vector<bool> cad_curve_list_selection;

    QLineEdit *cad_input, *arr_input;
    QTextEdit *oc_input;
    QComboBox *oc_method_box;
    Graphic_layer *axis;

    QWidget *central_widget;
    CGAL::Qt_widget *widget;
    CGAL::Qt_widget_standard_toolbar *stoolbar;

    int cur_method;
};

#endif


