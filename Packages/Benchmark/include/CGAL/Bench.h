// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>
#ifndef CGAL_BENCH_H
#define CGAL_BENCH_H

#include <time.h>
#include <signal.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip>

#if !(defined _MSC_VER)
#include <unistd.h>
#endif

CGAL_BEGIN_NAMESPACE

/*!
 */
class Bench_base {
public:
  /*!
   */
  static int get_name_length() { return m_name_length; }
  static void set_name_length(int nameLength) { m_name_length = nameLength; }
    
  /*!
   */
  static void print_header()
  {
    std::cout << std::setw(m_name_length) << std::left << "Bench"
              << " Bench    Ops      Total    Single   Num Ops \n";

    std::cout << std::setw(m_name_length) << std::left << "    Name"
              << "     Time      Num Ops Time  Op Time  Per Sec\n";

    std::cout << std::setw(m_name_length) << std::setfill('-') << ""
              << " -------- -------- -------- -------- --------\n"
              << std::setfill(' ');
  }
    
protected:
#if (!defined _MSC_VER)
  /*!
   */
  static void process_signal(int sig) { m_got_signal = true; }
#endif

  static bool m_got_signal;
  static int m_name_length;    
};

template <class Benchable>
class Bench : public Bench_base {
private:
  typedef Bench<Benchable> Self;
    
public:
  /*!
   */
  Bench(std::string name = "", int seconds = 1, bool print_header = true) :
      m_name(name), 
      m_print_header(print_header),
      m_header_printed(false),
      m_factor(1.0f), m_seconds(seconds),
      m_samples(0), m_period(0.0f), m_iterations(0)
  {}

  virtual ~Bench() {}
  
  int get_iterations() const { return m_iterations; }
  int get_seconds() const { return m_seconds; }
  int get_samples() const { return m_samples; }

  void set_iterations(int iterations) { m_iterations = iterations; }
  void set_seconds(int seconds) { m_seconds = seconds; }
  void set_samples(int samples) { m_samples = samples; }

  float get_period() const { return m_period; }

  /*!
   */
  void operator()()
  {
    if (m_benchable.init() < 0) return;
    loop();
    m_benchable.clean();
    print_results();
  }

  /*!
   */
  void loop()
  {
    int i;
    time_t cycle_secs = UINT_MAX / CLOCKS_PER_SEC;

    m_benchable.sync();
    if (m_samples == 0) {
      if (m_iterations != 0) {
        time_t time_start_secs = time(0);
        clock_t time_start_ticks = clock();
        for (i = 0; i < m_iterations; i++) m_benchable.op();
        m_benchable.sync();
        clock_t time_end_ticks = clock();
        time_t sec_end_secs = time(0);

        clock_t time_interval_ticks = time_end_ticks - time_start_ticks;
        time_t time_interval_secs = time_interval_ticks / CLOCKS_PER_SEC;
        time_t time_interval_approx_secs = sec_end_secs - time_start_secs;

        unsigned int num_cycles = time_interval_approx_secs / cycle_secs;

        time_t time_interval_tmp = num_cycles * cycle_secs + time_interval_secs;
        if (time_interval_tmp > time_interval_approx_secs) {
          time_t diff = time_interval_tmp - time_interval_approx_secs;
          if (diff > cycle_secs / 2)
            num_cycles--;
        } else {
          time_t diff = time_interval_approx_secs - time_interval_tmp;
          if (diff > cycle_secs / 2)
            num_cycles++;
        }
        
        float time_interval = num_cycles * cycle_secs +
          (float) time_interval_ticks / (float) CLOCKS_PER_SEC;
        m_samples = (int) ((float) m_seconds / time_interval);
        if (m_samples == 0) m_samples = 1;
        m_samples *= m_iterations;
      } else {
#if (defined _MSC_VER)
        std::cerr << "signal() is not supported by windows!" << std::endl
                  << "Set the number of iterations directly "
                  << "or set the number of samples." << std::endl;
#else
        m_got_signal = false;
        
        signal(SIGALRM, &Self::process_signal);
        alarm(m_seconds);
        do {
          m_benchable.op();
          m_samples++;
        } while (!m_got_signal);
#endif
      }
    } else {
      m_seconds = 0;
    }

    m_benchable.sync();
    time_t time_start_secs = time(0);
    clock_t time_start_ticks = clock();
    for (i = 0; i < m_samples; i++) m_benchable.op();
    m_benchable.sync();
    clock_t time_end_ticks = clock();
    time_t sec_end_secs = time(0);

    clock_t time_interval_ticks = time_end_ticks - time_start_ticks;
    time_t time_interval_secs = time_interval_ticks / CLOCKS_PER_SEC;
    time_t time_interval_approx_secs = sec_end_secs - time_start_secs;
    unsigned int num_cycles = time_interval_approx_secs / cycle_secs;
    time_t time_interval_tmp = num_cycles * cycle_secs + time_interval_secs;
    if (time_interval_tmp > time_interval_approx_secs) {
      time_t diff = time_interval_tmp - time_interval_approx_secs;
      if (diff > cycle_secs / 2)
        num_cycles--;
    } else {
      time_t diff = time_interval_approx_secs - time_interval_tmp;
      if (diff > cycle_secs / 2)
        num_cycles++;
    }
    m_period = num_cycles * cycle_secs +
      (float) time_interval_ticks / (float) CLOCKS_PER_SEC;

    // std::cout << m_samples << " iterations, " << m_period <<" seconds"
    //           << std::endl;
  }

  /*!
   */
  virtual void print_results()
  {
    fflush(stdout);
    if (!m_header_printed && m_print_header) print_header();
    m_header_printed = true;
    int count = (int) ((float) m_samples * m_factor);
    std::cout << std::setw(m_name_length) << std::left
              << m_name.substr(0, m_name_length).c_str() << " "
              << std::setw(8) << std::right << m_seconds << " "
              << std::setw(8) << std::right << count << " "
              << std::setw(8) << std::right << std::setprecision(4)
              << std::fixed << m_period << " "
              << std::setw(8) << std::right << std::setprecision(4)
              << std::fixed << m_period / (float) count << " "
              << std::setw(8) << std::right << std::setprecision(4)
              << std::fixed << (float) count / m_period
              << std::endl;
  }

  Benchable & get_benchable() { return m_benchable; }
    
private:
  Benchable m_benchable;

  std::string m_name;
  bool m_print_header;
  bool m_header_printed;
  float m_factor;
  int m_seconds;
  int m_samples;
  float m_period;
  int m_iterations;
};

CGAL_END_NAMESPACE

#endif
