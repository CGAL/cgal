// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_BENCHMARK_BENCHMARK_HPP
#define CGAL_BENCHMARK_BENCHMARK_HPP

/*! \file
 * A utility template class(es) to perform benchmarks.
 */

#include <time.h>
#include <signal.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <limits.h>

#if !(defined _MSC_VER)
#include <unistd.h>
#endif

#include "CGAL/Benchmark/config.hpp"

CGAL_BENCHMARK_BEGIN_NAMESPACE

/*! The non-template base class.
 * Contains the definition of static global data that 1. controls the
 * visual characteristics of the displayed fields, 2. displays the header,
 * and 3. handles signals.
 */
class Benchmark_base {
public:
  /*! Obtain the length of the name field */
  static int get_name_length() { return m_name_length; }

  /*! Set the the length of the name field */
  static void set_name_length(int nameLength) { m_name_length = nameLength; }
    
  /*! Display the header */
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
  /*! Mark the global flag that indicates that the alarm signal was recieved */
  static void process_signal(int sig) { m_got_signal = true; }
#endif

  /*! Indicates that alarm signal was recieved */
  static bool m_got_signal;

  /*! The length of the displayed name-field */
  static int m_name_length;    
};

/*! The main class.
 */
template <class Benchable>
class Benchmark : public Benchmark_base {
private:
  typedef Benchmark<Benchable> Self;
    
public:
  /*! Construcor */
  Benchmark(std::string name = "", int seconds = 1, bool print_header = true) :
      m_name(name), 
      m_print_header(print_header),
      m_header_printed(false),
      m_factor(1.0f), m_seconds(seconds),
      m_samples(0), m_period(0.0f), m_iterations(0)
  {}

  /*! Destructor */
  virtual ~Benchmark() {}
  
  int get_iterations() const { return m_iterations; }
  int get_seconds() const { return m_seconds; }
  int get_samples() const { return m_samples; }

  void set_iterations(int iterations) { m_iterations = iterations; }
  void set_seconds(int seconds) { m_seconds = seconds; }
  void set_samples(int samples) { m_samples = samples; }

  float get_period() const { return m_period; }

  /*! Perform the benchmark */
  void operator()()
  {
    if (m_benchable.init() < 0) return;
    loop();
    m_benchable.clean();
    print_results();
  }

  /*! Loop over the operation, the time of which is measured */
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

  /*! Print the results */
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

  /*! Obtain the benchable instance */
  Benchable & get_benchable() { return m_benchable; }
    
private:
  /*! The benchable instance */
  Benchable m_benchable;

  /*! The title */
  std::string m_name;

  /*! Indicates whether to print the header */
  bool m_print_header;

  /*! Indicates whether the header was printed already */
  bool m_header_printed;

  /*! Scales the number of operations (in case the measured operation already
   * consists of a sequence of many operations)
   */
  float m_factor;

  /*! The user-specified time in seconds the benchmark should approximately
   * last
   */
  int m_seconds;

  /*! The exact number of samples the operation being measured should be
   * executed based on the m_seconds estimate
   */
  int m_samples;

  /*! The time in seconds the loop over the operation being measured lasted */
  float m_period;

  /*! The iser-specifed number of times the operation being measured should be
   * executed. If this is set m_seconds has no effect.
   */
  int m_iterations;
};

CGAL_BENCHMARK_END_NAMESPACE

#endif
