#ifndef DIAL_OPT_H_
#define DIAL_OPT_H_

#include "ui_options.h"

class Dialog_options : public QDialog, private Ui::Dialog_options
{
    Q_OBJECT

public Q_SLOTS:
    void on_multiple_choice_checkbox_toggled(bool checked)
    {
      if (checked && get_mchoice() == 0)
        set_random_sample_size(10);
      else if (checked == false)
        set_random_sample_size(0);
    }

    void on_mchoice_spinbox_valueChanged(int new_v)
    {
      set_multiple_choice_checkbox(new_v != 0);
    }

public:
    Dialog_options(QWidget * = 0)
    {
        setupUi(this);
    }

    void set_all_ranges()
    {
        verbose_spinbox->setRange(0, 2);
        mchoice_spinbox->setRange(0, 500);
        percent_spinbox->setRange(0, 100);
        relocation_spinbox->setRange(0, 50);
        relevance_spinbox->setRange(0., 100000.);
    }

    int get_verbose() const { return verbose_spinbox->value(); }
    void set_verbose(int verbose) { verbose_spinbox->setValue(verbose); }

    int get_mchoice() const { return mchoice_spinbox->value(); }
    void set_random_sample_size(const int mchoice) { mchoice_spinbox->setValue(mchoice); }

    double get_percent() const { return percent_spinbox->value(); }
    void set_percent(const double percent) { percent_spinbox->setValue(percent); }

    int get_relocation() const { return relocation_spinbox->value(); }
    void set_relocation(const int value) {
      return relocation_spinbox->setValue(value);
    }

    double get_relevance() const { return relevance_spinbox->value(); }
    void set_relevance(double ghost) { relevance_spinbox->setValue(ghost); }

    bool get_use_flip() const { return use_flip_checkbox->isChecked(); }
    void set_use_flip(const bool flip) {
      return use_flip_checkbox->setChecked(flip);
    }

    bool get_multiple_choice_checkbox() const {
      return multiple_choice_checkbox->isChecked(); }
    void set_multiple_choice_checkbox(const bool mc) {
      return multiple_choice_checkbox->setChecked(mc);
    }

    double get_line_thickness() const { return thickness_spinbox->value(); }
    void set_line_thickness(const double t) { thickness_spinbox->setValue(t); }

    double get_point_size() const { return point_size_spinbox->value(); }
    void set_point_size(const double t) { point_size_spinbox->setValue(t); }

    double get_vertex_size() const { return vertex_size_spinbox->value(); }
    void set_vertex_size(const double t) { vertex_size_spinbox->setValue(t); }
};

#endif
