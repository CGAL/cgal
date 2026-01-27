#ifndef QT_PROGRESS_BAR_CALLBACK_H
#define QT_PROGRESS_BAR_CALLBACK_H

#include <CGAL/Real_timer.h>

#include <QProgressDialog>
#include <QApplication>

class Qt_progress_bar_callback
{

private:

  CGAL::Real_timer timer;
  double t_start;
  mutable double t_latest;

  mutable std::size_t nb;

  QProgressDialog* dialog;

public:

  Qt_progress_bar_callback()
  {
  }

  Qt_progress_bar_callback(const char* title, QWidget* parent)
    : dialog (new QProgressDialog (QString(title),
                                   QString("Cancel"),
                                   0, 100,
                                   parent))
  {
    dialog->setMinimumDuration(0);
    dialog->setWindowModality(Qt::WindowModal);
    timer.start();
    t_start = timer.time();
    t_latest = t_start;
  }
  ~Qt_progress_bar_callback()
  {
  }

  bool operator() (double advancement) const
  {
    // Avoid calling time() at every single iteration, which could
    // impact performances very badly
    ++ nb;
    if (advancement != 1 && nb % 1000 != 0)
      return true;

    // If the limit is reach, interrupt the algorithm
    double t = timer.time();
    if (advancement == 1 || (t - t_latest) > 0.25)
    {
      dialog->setValue (int (100. * advancement));
      t_latest = t;

      if (dialog->wasCanceled())
        return false;
    }

    return true;
  }
};


#endif // QT_PROGRESS_BAR_CALLBACK_H
