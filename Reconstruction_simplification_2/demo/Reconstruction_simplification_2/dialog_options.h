#ifndef DIAL_OPT_H_
#define DIAL_OPT_H_

#include "ui_options.h"

class Dialog_options : public QDialog, private Ui::Dialog_options
{
    Q_OBJECT

public:
    Dialog_options(QWidget *parent = 0)
    {
        setupUi(this);
    }
    
    void set_all_ranges()
    {
        verbose_spinbox->setRange(0, 2);
        mchoice_spinbox->setRange(0, 500);
        percent_spinbox->setRange(0, 100);
        norm_tol_spinbox->setRange(0., 500.);
        tang_tol_spinbox->setRange(0., 500.);
        alpha_spinbox->setRange(0., 100.);
        relocation_spinbox->setRange(0, 50);
        ghost_spinbox->setRange(0., 100000.);
    }
    
    int get_verbose() const { return verbose_spinbox->value(); }
    void set_verbose(int verbose) { verbose_spinbox->setValue(verbose); }    
    
    int get_mchoice() const { return mchoice_spinbox->value(); }
    void set_mchoice(const int mchoice) { mchoice_spinbox->setValue(mchoice); }
    
    double get_percent() const { return percent_spinbox->value(); }
    void set_percent(const double percent) { percent_spinbox->setValue(percent); }
    
    double get_norm_tol() const { return norm_tol_spinbox->value(); }
    void set_norm_tol(const double tol) { norm_tol_spinbox->setValue(tol); }
    
    double get_tang_tol() const { return tang_tol_spinbox->value(); }
    void set_tang_tol(const double tol) { tang_tol_spinbox->setValue(tol); }    

    double get_alpha() const { return alpha_spinbox->value(); }
    void set_alpha(const double alpha) { alpha_spinbox->setValue(alpha); }
    
    int get_relocation() const { return relocation_spinbox->value(); }
    void set_relocation(const int value) {
    	return relocation_spinbox->setValue(value);
    }
    
    double get_ghost() const { return ghost_spinbox->value(); }
    void set_ghost(double ghost) { ghost_spinbox->setValue(ghost); }
    
    bool get_use_flip() const { return use_flip_checkbox->isChecked(); }
    void set_use_flip(const bool flip) {
    	return use_flip_checkbox->setChecked(flip);
    }
    
    double get_line_thickness() const { return thickness_spinbox->value(); }
    void set_line_thickness(const double t) { thickness_spinbox->setValue(t); }

    double get_point_size() const { return point_size_spinbox->value(); }
    void set_point_size(const double t) { point_size_spinbox->setValue(t); }

    double get_vertex_size() const { return vertex_size_spinbox->value(); }
    void set_vertex_size(const double t) { vertex_size_spinbox->setValue(t); }    
};

#endif
