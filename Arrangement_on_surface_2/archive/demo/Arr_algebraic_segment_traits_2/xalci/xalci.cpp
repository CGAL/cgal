// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/).
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
// ----------------------------------------------------------------------------
//
// Library       : AlciX
// File          : demos/xalci/xalci.cpp
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.43 $
// Revision_date : $Date: 2009-06-30 13:14:58 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-inf.mpg.de>
//                 Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#define NDEBUG 1

//#define CGAL_NO_LEDA

#include <CGAL/basic.h>
#include <CGAL/Timer.h>

#include <set>
#include <ctime>
#include <fstream>

#define CGAL_POLYNOMIAL_USE_NTL_MUL

#include "include/misc.h"
#include "include/xalci.h"

#define XALCI_USE_FLAT_COLOR_SCHEME 0

#ifdef CGAL_ACK_BENCHMARK_RES
namespace CGAL {
CGAL::Timer res_tm;
}
#endif

static const Kernel_2& kernel_2 = CKvA_2::instance().kernel();

bool check_testsuite(int argc, char** argv) {
    if (argc > 1) {
        for (int a = 0; a < argc; a++) {
            if (std::strcmp( argv[a], "--test-suite") == 0) {
                std::cerr << "This interactive program terminates "
                          << "immediately for the test-suite." << std::endl;
                return true;
            }
        }
    }
    return false;
}

#if !CGAL_HAS_QT3

int main (int argc, char** argv) {

    if(check_testsuite(argc, argv)) 
        return 0;
    
    std::cerr << "This demo requires Qt!" << std::endl;
    return 0;
}
#else // !CGAL_USE_QT

extern QColor rasterize_colors[];
extern int n_rast_colors;

extern Object_vector cad_objects, oc_objects, arr_objects;
extern Object_vector* curr_objects;
extern bool subdiv_layer_changed;
 
extern Graphic_layer *subdiv_layer;
extern Layers cad_layers, oc_layers, arr_layers;

extern CGAL::Timer timer;

typedef std::map<std::size_t, int> Color_map;

// picks up a color index based on id parameter
static int pick_color(std::size_t id, Color_map& color_map) {

    int cindex = color_map.size();
    std::pair<Color_map::iterator, bool> ret =
        color_map.insert(std::make_pair(id, cindex));
    if(!ret.second) // already exists in the map
        cindex = ret.first->second;
    return cindex;
}

void xAlci_main_window::cad_curve_list_click() {

    for(int i = 0; i < static_cast<int>(cad_curve_list->count()); i++) {
        if(cad_curve_list->isSelected(i) && !cad_curve_list_selection[i]) {

            cad_curve_list_selection[i] = true;
            for(int j=0;j<static_cast<int>(cad_seg_list->count());j++) {
                int fidx = cad_layers[j]->get_first_index(),
                    sidx = cad_layers[j]->get_second_index();

                if(fidx == i || sidx == i)
                    if(cad_curve_list_selection[fidx] &&
                            cad_curve_list_selection[sidx]) 
                        cad_seg_list->setSelected(j, true);
            }
        }
        if(!cad_curve_list->isSelected(i) && cad_curve_list_selection[i]) {
                cad_curve_list_selection[i]=false;
            for(int j = 0; j < static_cast<int>(cad_seg_list->count()); j++) {
                if(cad_layers[j]->get_first_index() == i ||
                        cad_layers[j]->get_second_index() == i) 
                    cad_seg_list->setSelected(j,false);          
            }
        }
    }
}

void xAlci_main_window::oc_switch_method(int index)
{
    if(index != cur_method) {
        cur_method = index;
        if(cur_method == 0) { // segment renderer
            subdiv_layer->deactivate();
            oc_activate_layers();
            oc_analyse_btn->setEnabled(true);
            oc_seg_list->setEnabled(true);
            oc_complete_check->setEnabled(true);
        } else {            // space subdivision
            oc_deactivate_layers();
            subdiv_layer->activate();
            oc_analyse_btn->setEnabled(false);
                        oc_seg_list->setEnabled(false);
            oc_complete_check->setEnabled(false);
            oc_rasterize_btn->setEnabled(true);
        }
        widget->redraw();
    }
}

void xAlci_main_window::arr_analyse_click() {

  arr_rasterize_btn->setEnabled(false);
  
  //CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
  //std::scientific(std::cout);
  //std::cout.precision(10);
  
  std::vector<Poly_int2> temp_arr_curves;
  if(read_polys_from_file(arr_input->text(),
                          std::back_inserter(temp_arr_curves))) {
    if(arr_method->selectedId()==-1) {
      QMessageBox::warning(this,"Error","No Sweep-line method specified",0);
      return ;
    }
    arr_curves = temp_arr_curves;
    arr_compute_arrangement();
  }
  if(arr_complete_check->isChecked()) {
    arr_rasterize_btn->setEnabled(true);
  }
}

void xAlci_main_window::oc_analyse_click()
{
    int i;
    oc_seg_list->clear();
    oc_objects.clear();
    for(Layers::iterator it= oc_layers.begin(); it != oc_layers.end(); it++) {
        widget->detach(*it);
        delete (*it);
    }
    
    oc_layers.clear();
    oc_rasterize_btn->setEnabled(false);
    
    Poly_int2 f;
    if(!input_poly(f, oc_input->text().ascii()))
        return;

//         *rational_fx = typename RP_traits::Differentiate()(*rational_y, 0);
//     *rational_fy = typename RP_traits::Differentiate()(*rational_y, 1);

    CGAL::Polynomial_traits_d< Poly_int2 >::Differentiate diff;

    Poly_int2 fx = diff(f, 0), fxx = diff(f, 0), fxy = diff(fx, 1), 
        fy = diff(f, 1), fyy = diff(fy, 1);
    Poly_int2 ress = fxx*fy*fy - ((fx*fy*fxy)*Poly_int1(2,0)) + fyy*fx*fx,
        res1 = f*fx, res2 = f*fy;
  
    CGAL::set_pretty_mode(std::cout);
    std::cout << "curv:\n " << ress << "\n\n";
    std::cout << "fx:\n " << fx << "\n\n";
    std::cout << "fy:\n " << fy << "\n\n";
    std::cout << "f*fx:\n " << res1 << "\n\n";
    std::cout << "f*fy:\n " << res2 << "\n\n";

    CGAL::set_ascii_mode(std::cout);
    std::cout << "f:\n " << f << "\n\n";
    CGAL::set_pretty_mode(std::cout);
    std::cout << "f:\n " << f << "\n\n";
    
    timer.reset();
    timer.start();
    
//     CGAL::res_tm.reset();

    {

    Kernel_2::Construct_curve_2 cc_2 = kernel_2.construct_curve_2_object();
    
Lbegin:
    try {
        Curve_analysis_2 curve = cc_2(f);
        CKvA_2::instance().make_x_monotone_2_object()(curve,
            std::back_inserter(oc_objects));

    } catch(...) {
        std::cout << "Seemingly non-square free.. restarting..\n";
        f = make_square_free(f);

        std::cout << "Square-free part: " << f << "\n";
        goto Lbegin;
    }

        timer.stop();
        std::cout << "\n\nAnalyse elapsed time: " << timer.time() << std::endl;
//         std::cout << "\nResultant time: " << CGAL::res_tm.time() << "\n";
        std::cout << oc_objects.size() << " arcs found (incl isolated points)"
             << std::endl;
        
        Object_vector::const_iterator oit;
        Arc_2 arc;
        Point_2 pt;
        
        for(oit = oc_objects.begin(), i = 0; oit != oc_objects.end();
                oit++, i++) {

            std::ostringstream os;
            if(CGAL::assign(arc, *oit))
                print_arc(arc, os);
            else if(CGAL::assign(pt, *oit)) {
                print_point(pt, os);
                os << "; isolated";
            } else
                CGAL_error_msg("Malformed object found..\n");

            oc_seg_list->insertItem(os.str());
#if !XALCI_USE_FLAT_COLOR_SCHEME
            oc_layers.push_back(new Graphic_layer(widget, i, i));
#else
            oc_layers.push_back(new Graphic_layer(widget, i, i));
#endif
            oc_layers[i]->deactivate();
        }
        //subdiv_renderer.set_polynomial(f);
    
    } 

    subdiv_layer_changed = true;
    //layers.push_back(new Graphic_layer(-1,widget));
    oc_rasterize_btn->setEnabled(oc_complete_check->isChecked());
    widget->redraw();
}

void xAlci_main_window::cad_partial_selection_click()
{
    std::vector<Poly_int2> temp_cad_curves;
    if(read_polys_from_file(cad_input->text(),
        std::back_inserter(temp_cad_curves))) {
        
        cad_curves = temp_cad_curves;
        std::stringstream strstr;
        strstr << "Choose polynomials (" << cad_curves.size() 
               << " in total)";
        curve_selection_dialog = new
             Curve_selection_dialog(central_widget,tr(strstr.str().c_str()),
                false, cad_curves.begin(), cad_curves.end());
            
        curve_selection_dialog->setCaption(tr(strstr.str().c_str()));
        int sel_res = curve_selection_dialog->exec();
    
        if(sel_res==QDialog::Accepted) {
    
            QListBox* curve_selection = curve_selection_dialog->curve_list;
            std::vector<Poly_int2>::iterator curve_it = cad_curves.begin();
            int n = static_cast<int>(cad_curves.size());
            CGAL_assertion(n==curve_selection->numRows());
        
            for(int i=0;i<n;i++) {
                if(!curve_selection->isSelected(i)) 
                    curve_it = cad_curves.erase(curve_it);
                else 
                    curve_it++;
            }
            CGAL_assertion(curve_it == cad_curves.end());
            cad_to_segments();
        }
    }
}


void xAlci_main_window::arr_partial_selection_click()
{
    std::vector<Poly_int2> temp_arr_curves;
  
    if(read_polys_from_file(arr_input->text(),
         std::back_inserter(temp_arr_curves))) {
        
        if(arr_method->selectedId() == -1) {
            QMessageBox::warning(this, "Error",
                "No Sweep-line method specified", 0);
            return;
        }

        arr_curves = temp_arr_curves;
        std::stringstream strstr;
        strstr << "Choose polynomials (" << arr_curves.size() << " in total)";
        
        curve_selection_dialog = new
            Curve_selection_dialog(central_widget, tr(strstr.str().c_str()),
                false, arr_curves.begin(), arr_curves.end());

        curve_selection_dialog->setCaption(tr(strstr.str().c_str()));

        int sel_res = curve_selection_dialog->exec();
        if(sel_res==QDialog::Accepted) {
            QListBox* curve_selection = curve_selection_dialog->curve_list;

            std::vector<Poly_int2>::iterator curve_it = arr_curves.begin();
            int n = static_cast<int>(arr_curves.size());
            CGAL_assertion(n==curve_selection->numRows());

            for(int i = 0; i < n; i++) {
                if(!curve_selection->isSelected(i)) 
                    curve_it = arr_curves.erase(curve_it);
                else 
                    curve_it++;
            }
            CGAL_assertion(curve_it == arr_curves.end());
            arr_compute_arrangement();
        }
    }
}


template<typename OutputIterator>
bool xAlci_main_window::read_polys_from_file(QString filename,
        OutputIterator out) {

    if(filename.isEmpty()) {
        QMessageBox::warning(this, "Error", "No input file specified", 0);
        return false;
    }

    std::ifstream ifstr(filename.ascii());
    if(!ifstr.is_open()) {
        std::stringstream err_msg;
        err_msg << "File \"" << filename.ascii() <<"\" not found" ;
        QMessageBox::warning(this, tr("Error"), tr(err_msg.str().c_str()), 0);
        return false;
    }

    Poly_int2 f;
    while(!ifstr.eof()) {
        int g = ifstr.get();
        if(g=='P') {
            ifstr.putback(g);
            ifstr >> f; 
            *out++ = f;
        }
    }
    return true;
}

void xAlci_main_window::cad_to_segments() {

    cad_seg_list->clear();
    cad_curve_list->clear();
    cad_objects.clear();

    for(Layers::iterator it = cad_layers.begin(); it != cad_layers.end();
            it++) {
        widget->detach(*it);
        delete (*it);
    } 
    cad_layers.clear();

    typedef std::vector<Curve_analysis_2> Curve_vector;
    Curve_vector curves;
    Object_vector arcs, intersections;

    std::vector<std::pair<int, int> > arc_info, inter_info, s_arcs_info;
    int curve_count = 0;
        
    typedef Kernel_2::Curve_pair_analysis_2 Curve_pair_analysis_2;
    typedef Curve_pair_analysis_2::Status_line_1 Status_line_1;

    Kernel_2::Construct_curve_2 cc_2 = kernel_2.construct_curve_2_object();
    Kernel_2::Construct_curve_pair_2 ccp_2 =
            kernel_2.construct_curve_pair_2_object();
    CKvA_2::Make_x_monotone_2 make_x_monotone(&CKvA_2::instance());
    
    int n = static_cast<int>(cad_curves.size());
    int number_of_steps = n + n*(n-1)/2;

    timer.reset();
    timer.start();
    
    QProgressDialog* progress = new QProgressDialog(
        tr("Please wait"), tr("Abort"), number_of_steps, this,
             "Progress", true);
    
    progress->setCaption("Analyse curves");
    connect(progress,SIGNAL(canceled),SLOT(cancel));
    progress->show();

    int progress_count=0;
    for(int j = 0; j < n; j++) {
        if(progress->wasCanceled())
                break;

        Curve_analysis_2 curve = cc_2(make_square_free(cad_curves[j]));
        progress->setProgress(progress_count++);
            
        int arc_before = arcs.size(), arc_no;
        make_x_monotone(curve, std::back_inserter(arcs));
        arc_no = arcs.size() - arc_before;

        std::stringstream sstr;
        sstr << "Curve no. " << curve_count << ", " << arc_no << "arcs";
        cad_curve_list->insertItem(sstr.str());
        cad_curve_list_selection.push_back(false);

        for(int i = 0; i < arc_no; i++)
            arc_info.push_back(std::make_pair(curve_count,curve_count));
        
//         int curve_int_count=0; 
        for(Curve_vector::const_iterator cit = curves.begin();
                cit != curves.end(); cit++) {

 //#warning "cad2segments is temporary (or permanently ?) out of service"
/*            if(progress->wasCanceled())
                break;
                
            Curve_pair_analysis_2 cpa;
            try {
                cpa = ccp_2(curve, *cit);
            } catch( ... ) {

                if(curve.polynomial_2() == cit->polynomial_2()) {
                    progress_count++;
                    continue;
                }
                
                Poly_int2 no_content, q, r, pair_gcd = NiX::gcd(
                    curve._internal_curve().f_primitive(),
                    cit->_internal_curve().f_primitive());
                    
                Poly_int1 d, pair_content = NiX::gcd(
                    curve._internal_curve().content(),
                    cit->_internal_curve().content());
                    
                Poly_int2::euclidean_division(cit->polynomial_2(),
                    pair_gcd, q, r);

                Curve_analysis_2 coprim_it = cc_2(q);
                cpa = ccp_2(curve, coprim_it);
            }

            for(int i = 0; i < cpa.number_of_status_lines_with_event(); i++) {

                Curve_pair_analysis_2::Status_line_1 cpa_line =
                    cpa.status_line_at_event(i);
                
                if(!cpa_line.is_intersection())
                    continue;

                X_coordinate_1 x0 = cpa_line.x();
                Kernel_2::Curve_analysis_2::Status_line_1 sline =
                    curve.status_line_for_x(x0);

                for(int j = 0; j < sline.number_of_events(); j++) {
                
                    int pair_idx = cpa_line.event_of_curve(j, 0);
                    Status_line_1::Arc_pair pair =
                        cpa_line.curves_at_event(pair_idx);

                    if(pair.first != -1 && pair.second != -1) {
                        // TODO: use Status_line_CA_1::algebraic_real_2() ?
                        Point_2 pt(x0, cpa.curve_analysis(0), j);
                        intersections.push_back(CGAL::make_object(pt));
                        inter_info.push_back(std::make_pair(
                            curve_count, curve_int_count));
                    }
                }
            }
            progress->setProgress(progress_count++);
            curve_int_count++;*/
        }
        curves.push_back(curve);
        curve_count++;
    }

    if(progress->wasCanceled()) {
        err_msg_dialog->setCaption("User information");
        err_msg_dialog->message("Computation has been aborted. \
             Results are still presented, but use at your own risk!");
    } else {
        CGAL_assertion(progress_count == number_of_steps);
    }
  
    delete progress;
    timer.stop();

    std::cout << "\n\nAnalyse elapsed time: " << timer.time() << std::endl;
    std::cout << arcs.size() << " arcs found" << std::endl;
    
    CGAL_assertion(arcs.size() == arc_info.size());
    
    std::cout << intersections.size() << " intersections found" << std::endl;
    CGAL_assertion(intersections.size() == inter_info.size());
    
    std::copy(arcs.begin(), arcs.end(), std::back_inserter(cad_objects));
    
    std::copy(arc_info.begin(), arc_info.end(),
        std::back_inserter(s_arcs_info));

    std::copy(intersections.begin(), intersections.end(),
        std::back_inserter(cad_objects));
        
    std::copy(inter_info.begin(), inter_info.end(),
              std::back_inserter(s_arcs_info));

    Color_map color_map;
        
    int i = 0;
    Point_2 pt;
    Arc_2 arc;
    
    for(Object_vector::const_iterator oit = cad_objects.begin();
            oit != cad_objects.end(); oit++, i++) {

        int cindex = 0;
        std::ostringstream os;
        if(CGAL::assign(arc, *oit)) {
            print_arc(arc, os);
            cindex = arc.curve().id();
        } else if(CGAL::assign(pt, *oit)) {
            print_point(pt, os);
            os << "; isolated";
        } else
            CGAL_error_msg("Malformed object..\n");

        cad_seg_list->insertItem(os.str());
        cindex = pick_color(cindex, color_map);
        cad_layers.push_back(new Graphic_layer(widget, i, cindex,
             s_arcs_info[i].first, s_arcs_info[i].second));
        cad_layers[i]->deactivate();
    }

   
    cad_rasterize_btn->setEnabled(cad_complete_check->isChecked());
    widget->redraw();
}

void xAlci_main_window::arr_compute_arrangement() {

    arr_node_list->clear();
    arr_edge_list->clear();
    arr_objects.clear();
    for(Layers::iterator it= arr_layers.begin(); it != arr_layers.end();
            it++) {
        widget->detach(*it);
        delete (*it);
    } 
    arr_layers.clear();
    
    timer.reset();
    timer.start();
    
    Kernel_2::Construct_curve_2 cc_2 = kernel_2.construct_curve_2_object();
    typedef std::vector<Curve_analysis_2> Curve_vector;

    Curve_vector curves;

    int n = static_cast<int>(arr_curves.size());
    for(int j = 0; j < n; j++) {

        Curve_analysis_2 tmp = cc_2(make_square_free(arr_curves[j]));
        //CKvA_2::instance().make_x_monotone_2_object()(tmp,
            //  std::back_inserter(oc_objects));
        curves.push_back(tmp);
    }

    if(arr_method->selectedId()==0) {
#if 0
#ifdef CGAL_USE_LEDA
    
      typedef Arrangement_2::SoX_Arrangement_2 Leda_graph_2;
            
      Leda_graph_2 leda_graph_2;
      
      arrangement(segments.begin(),
                  segments.end(),
                  leda_graph_2,
                  true,
                  GAPS_2_inst());
      
      timer.stop();
      
      std::cout << "\n\nAnalyse elapsed time: " << timer.time() << std::endl;
      int i=0;
      std::stringstream strstr1;
      strstr1 << "<b>Nodes (" << leda_graph_2.number_of_nodes() << " in total):</b>";
      arr_node_label->setText(strstr1.str());
      arr_node_label->update();
      std::stringstream strstr2;
      strstr2 << "<b>Edges (" << leda_graph_2.number_of_edges()/2 << " in total):</b>";
      arr_edge_label->setText(strstr2.str());
      arr_edge_label->update();
      Leda_graph_2::edge e;
      Leda_graph_2::node v;
      forall_nodes(v,leda_graph_2) {
        std::stringstream buf;
        Point_2 p = leda_graph_2.inf(v);
        buf << "At: ";
        out_point(p,buf);
        arr_node_list->insertItem(buf.str());
        arr_layers.push_back(new Graphic_layer(i,widget,static_cast<QObject*>(0),static_cast<const char*>(0),random_color()));
        arr_layers[i]->deactivate(); 
        arr_objects.push_back(Arc_2(p));
        i++;
      }
      forall_edges(e,leda_graph_2) {
        std::stringstream buf;
        Arc_2 too_long_seg = leda_graph_2.inf(e);

        Point_2 source = leda_graph_2.inf(leda_graph_2.source(e)),
          target = leda_graph_2.inf(leda_graph_2.target(e));

        //CGAL_assertion(source.curve().is_identical(too_long_seg.support()));
        //CGAL_assertion(target.curve().is_identical(too_long_seg.support()));
        Arc_2 seg = too_long_seg.trim(source,target);
        if(std::find(arr_objects.begin(),
                     arr_objects.end(),
                     seg)==arr_objects.end()) {
        
          //if(true) { 
          buf << "from: ";  
          out_point(seg.source(), buf);
          buf << " to: ";
          out_point(seg.target(), buf);
          if(!seg.is_vertical())
            buf << " segment arcno: " << seg.arcno();
          else
            buf << " vertical";
          arr_edge_list->insertItem(buf.str());
          arr_layers.push_back(new Graphic_layer(i,widget,static_cast<QObject*>(0),static_cast<const char*>(0),random_color()));
          arr_layers[i]->deactivate(); 
          arr_objects.push_back(seg);
          i++;
        }
      }
      
#endif // CGAL_USE_LEDA
#endif
    } else if(arr_method->selectedId()==1) {
#if 1
        typedef CGAL::Arrangement_2<CKvA_2> Arrangement;
        Arrangement arr;

        CGAL::insert(arr, curves.begin(), curves.end(), boost::false_type());
                  
        timer.stop();
        std::cout << "\n\nAnalyse elapsed time: " << timer.time() << std::endl;

        std::stringstream outs1, outs2;
        outs1 << "<b>Nodes (" << (arr.number_of_vertices() -
            arr.number_of_isolated_vertices()) << " in total):</b>";
        arr_node_label->setText(outs1.str());
        arr_node_label->update();

        outs2 << "<b>Edges (" << (arr.number_of_edges() +
            arr.number_of_isolated_vertices()) << " in total):</b>";
        arr_edge_label->setText(outs2.str());
        arr_edge_label->update();
        
        Color_map color_map;
        
        int i = 0, cindex;
        for(Arrangement::Vertex_const_iterator vit = arr.vertices_begin();
                vit != arr.vertices_end(); vit++) {

            std::ostringstream os;
            const Point_2& pt = vit->point();
            print_point(pt, os);

            if(!vit->is_isolated()) {
                arr_node_list->insertItem(os.str());
                continue;
            }
            
            os << "; isolated";
            arr_edge_list->insertItem(os.str());
#if !XALCI_USE_FLAT_COLOR_SCHEME
            cindex = i;            
#else
            cindex = pick_color(pt.curve().id(), color_map);
#endif            
            
            arr_layers.push_back(new Graphic_layer(widget, i, cindex));
            arr_layers[i]->deactivate();
            arr_objects.push_back(CGAL::make_object(pt));
            i++;
        }

        for(Arrangement::Edge_const_iterator eit = arr.edges_begin();
                eit != arr.edges_end(); eit++, i++) {

            std::ostringstream os;
            const Arc_2& arc = eit->curve();
            print_arc(arc, os);
            arr_edge_list->insertItem(os.str());

#if !XALCI_USE_FLAT_COLOR_SCHEME            
            cindex = i;
#else
            cindex = pick_color(arc.curve().id(), color_map);
#endif

            arr_layers.push_back(new Graphic_layer(widget, i, cindex));
            arr_layers[i]->deactivate();
            arr_objects.push_back(CGAL::make_object(arc));
       }
#endif
   }

    arr_rasterize_btn->setEnabled(cad_complete_check->isChecked());
    widget->redraw();
}

Poly_int2 xAlci_main_window::make_square_free(const Poly_int2& poly) {
#if !AcX_SQRT_EXTENSION

    return kernel_2.decompose_2_object()(poly);
#else
    return poly; // no gcds available for sqrt_extensions
#endif
}

void xAlci_main_window::visualize()
{
    show();
    widget->zoom(1);
}


void xAlci_main_window::cad_analyse_click()
{
    cad_rasterize_btn->setEnabled(false);
  
//     CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
    //std::scientific(std::cout);
    //std::cout.precision(10);
    std::vector<Poly_int2> temp_cad_curves;
    if(read_polys_from_file(cad_input->text(),
            std::back_inserter(temp_cad_curves))) {

            
        cad_curves = temp_cad_curves;
        cad_to_segments();
    }
}

void xAlci_main_window::tab_changed(QWidget*) {

    Layers *detach[2], *attach;
    
    if(tab_widget->currentPage() == cad_tab) {
        curr_objects = &cad_objects;
        attach = &cad_layers;
        detach[0] = &oc_layers;
        detach[1] = &arr_layers;
        
    } else if(tab_widget->currentPage() == one_curve_tab) {
        curr_objects = &oc_objects;
        attach = &oc_layers;
        detach[0] = &cad_layers;
        detach[1] = &arr_layers;
        
    } else {
        curr_objects = &arr_objects;
        attach = &arr_layers;
        detach[0] = &cad_layers;
        detach[1] = &oc_layers;
    }

    Layers::iterator lit;
    for(int i = 0; i < 2; i++)
        for(lit = detach[i]->begin(); lit != detach[i]->end(); lit++) {
            widget->detach(*lit);
        }

    for(lit = attach->begin(); lit != attach->end(); lit++)
        widget->attach(*lit);
}

#include "xalci.moc"

int main(int argc, char **argv)
{
    if (check_testsuite(argc, argv)) {
        return (0);
    }
    QApplication app(argc, argv);
    xAlci_main_window widget(1024, 700); // physical window size
    app.setMainWidget(&widget);
    widget.setCaption("Curve rasterizer demo");
    widget.setMouseTracking(TRUE);
    widget.visualize();
    return app.exec();
}

#endif // CGAL_HAS_QT3
