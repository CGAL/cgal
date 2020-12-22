#ifndef PROGRESS_BAR_CALLBACK_H
#define PROGRESS_BAR_CALLBACK_H

#include <CGAL/Real_timer.h>

struct Progress_bar_callback
{
  mutable std::size_t nb;
  CGAL::Real_timer timer;
  double t_start;
  mutable double t_start_estimate;
  mutable double t_latest;
  int bar_size;
  mutable std::size_t string_size;

  Progress_bar_callback(const char* title = NULL,
                        void* = NULL) // Hack to have same constructor as Qt_progress_bar_callback
    : nb(0), bar_size (30), string_size(0)
  {
    if (title != NULL)
      std::cerr << title << std::endl;

    timer.start();
    t_start = timer.time();
    t_start_estimate = 0.;
    t_latest = t_start;
  }

  bool operator()(double advancement) const
  {
    // Avoid calling time() at every single iteration, which could
    // impact performances very badly
    ++ nb;
    if (advancement != 1 && nb % 1000 != 0)
      return true;

    // If the limit is reach, interrupt the algorithm
    double t = timer.time();
    if (advancement == 1 || (t - t_latest) > 1.)
    {
      if (advancement != 0. && t_start_estimate == 0.)
        t_start_estimate = t_latest;

      std::ostringstream oss;
      oss << "[";
      int adv = int(advancement * bar_size);
      for (int i = 0; i < adv; ++ i)
        oss << "=";
      if (adv != bar_size)
        oss << ">";
      for (int i = adv; i < bar_size; ++ i)
        oss << " ";

      oss << "] " << int(advancement * 100) << "% (";

      display_time (oss, t);

      oss << " elapsed, ";

      display_time (oss, estimate_remaining(t, advancement));

      oss << " remaining)";

      std::string bar_string = oss.str();
      std::cerr << "\r" << bar_string;
      for (std::size_t i = bar_string.size(); i < string_size; ++ i)
        std::cerr << " ";
      string_size = (std::max) (string_size, bar_string.size());

      t_latest = t;

      if (advancement == 1)
        std::cerr << std::endl;
    }

    return true;
  }

  void display_time (std::ostringstream& oss, double seconds) const
  {
    if (seconds == std::numeric_limits<double>::infinity())
    {
      oss << "unknown";
      return;
    }
    if (seconds > 3600.)
      {
        int hours = int(seconds / 3600.);
        oss << hours << "h";
        seconds -= hours * 3600.;
      }
    if (seconds > 60.)
      {
        int minutes = (int)(seconds / 60.);
        oss << minutes << "min";
        seconds -= minutes * 60.;
      }
    oss << int(seconds) << "sec";
  }

  double estimate_remaining (double seconds, double advancement) const
  {
    if (advancement == 0.)
      return std::numeric_limits<double>::infinity();

    return ((1. - advancement) * (seconds - t_start_estimate) / advancement);
  }
};


#endif // PROGRESS_BAR_CALLBACK_H
