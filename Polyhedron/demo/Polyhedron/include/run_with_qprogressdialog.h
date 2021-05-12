#ifndef RUN_WITH_QPROGRESSDIALOG_H
#define RUN_WITH_QPROGRESSDIALOG_H

#include <QProgressDialog>
#include <CGAL/Real_timer.h>
#include <CGAL/thread.h>

#include "Callback_signaler.h"

typedef CGAL::Parallel_if_available_tag Concurrency_tag;

class Signal_callback
{
private:

  CGAL::Real_timer timer;
  double t_start;
  mutable double t_latest;
  mutable std::size_t nb;

public:
  boost::shared_ptr<double> latest_adv;
  boost::shared_ptr<bool> state;
  boost::shared_ptr<Callback_signaler> signaler;

  Signal_callback(bool)
    : latest_adv (new double(0))
    , state (new bool(true))
    , signaler (new Callback_signaler())
  {
    timer.start();
    t_start = timer.time();
    t_latest = t_start;
  }
  ~Signal_callback()
  {
  }

  bool operator() (double advancement) const
  {
    if (!state)
      return false;

    *latest_adv = advancement;

    // Avoid calling time() at every single iteration, which could
    // impact performances very badly
    ++ nb;
    if (advancement != 1 && nb % 1000 != 0)
      return *state;

    // If the limit is reach, interrupt the algorithm
    double t = timer.time();
    if (advancement == 1 || (t - t_latest) > 0.25)
    {
      signaler->emit_signal (int (100. * advancement));

      t_latest = t;

      if (signaler->is_canceled)
        *state = false;
    }

    return *state;
  }

};


class Functor_with_signal_callback
{
protected:
  boost::shared_ptr<Signal_callback> m_callback;
public:

  Signal_callback* callback() { return m_callback.get(); }

  Functor_with_signal_callback()
    : m_callback (new Signal_callback(true)) { }

  virtual void operator()() = 0;
};


template <typename Functor>
void run_with_qprogressdialog (Functor& functor,
                               const char* title,
                               QWidget* mainWindow)
{
  return run_with_qprogressdialog<Concurrency_tag> (functor, title, mainWindow);
}

template <typename ConcurrencyTag, typename Functor>
void run_with_qprogressdialog (Functor& functor,
                               const char* title,
                               QWidget* mainWindow)
{
  mainWindow->setEnabled(false);
  QProgressDialog progress (QString(title),
                            QString("Cancel"),
                            0, 100,
                            mainWindow);
  progress.setMinimumDuration(0);

  Signal_callback* signal_callback = functor.callback();

  QEventLoop::connect (signal_callback->signaler.get(), SIGNAL(progressChanged(int)),
                       &progress, SLOT(setValue(int)));
  QEventLoop::connect (&progress, SIGNAL(canceled()),
           signal_callback->signaler.get(), SLOT(cancel()));

#ifdef CGAL_HAS_STD_THREADS
  if (boost::is_convertible<ConcurrencyTag, CGAL::Parallel_tag>::value)
  {
    CGAL::cpp11::thread thread (functor);

    while (*signal_callback->latest_adv != 1. &&
           *signal_callback->state)
    {
      CGAL::cpp11::sleep_for (0.1);
      QApplication::processEvents ();
    }

    thread.join();
  }
  else
#endif // Sequential version
  {
    progress.setWindowModality(Qt::WindowModal);
    functor();
  }

  mainWindow->setEnabled(true);
}


#endif //  RUN_WITH_QPROGRESSDIALOG_H
