#ifndef CGAL_BENCH_H
#define CGAL_BENCH_H

#include <time.h>
#include <signal.h>
#include <stdio.h>
#include <string>

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
  static int getNameLength() { return m_nameLength; }
  static void setNameLength(int nameLength) { m_nameLength = nameLength; }
    
  /*!
   */
  static void printHeader()
  {
    std::cout << std::setw(m_nameLength) << std::left << "Bench"
              << " Bench    Ops      Total    Single   Num Ops \n";

    std::cout << std::setw(m_nameLength) << std::left << "    Name"
              << "     Time      Num Ops Time  Op Time  Per Sec\n";

    std::cout << std::setw(m_nameLength) << std::setfill('-') << ""
              << " -------- -------- -------- -------- --------\n"
              << std::setfill(' ');
  }
    
protected:
#if (!defined _MSC_VER)
  /*!
   */
  static void processSignal(int sig) { m_gotSignal = true; }
#endif

  static bool m_gotSignal;
  static int m_nameLength;    
};

template <class Bench_user>
class Bench : public Bench_base {
private:
  typedef Bench<Bench_user> Self;
    
public:
  /*!
   */
  Bench(std::string name = "", int seconds = 1, bool printHeader = true) :
      m_name(name), 
      m_printHeader(printHeader),
      m_headerPrinted(false),
      m_factor(1.0f), m_seconds(seconds),
      m_samples(0), m_period(0.0f), m_iterations(0)
  {}
    
  int getIterations() const { return m_iterations; }
  int getSeconds() const { return m_seconds; }
  int getSamples() const { return m_samples; }

  void setIterations(int iterations) { m_iterations = iterations; }
  void setSeconds(int seconds) { m_seconds = seconds; }
  void setSamples(int samples) { m_samples = samples; }

  /*!
   */
  void operator()()
  {
    if (m_user.init() < 0) return;
    loop();
    m_user.clean();
    printResults();
  }

  /*!
   */
  void loop()
  {
    int i;

    m_user.sync();
    if (m_samples == 0) {
      if (m_iterations != 0) {
        clock_t time_start = clock();
        for (i = 0; i < m_iterations; i++) m_user.op();
        m_user.sync();
        clock_t time_end = clock();
        float time_interval = (float) (time_end - time_start) /
          (float) CLOCKS_PER_SEC;
        m_samples = (int) ((float) m_seconds / time_interval);
        if (m_samples == 0) m_samples = 1;
        m_samples *= m_iterations;
      } else {
#if (defined _MSC_VER)
        std::cerr << "signal() is not supported by windows!" << std::endl
                  << "Set the number of iterations directly "
                  << "or set the number of samples." << std::endl;
#else
        m_gotSignal = false;
        
        signal(SIGALRM, &Self::processSignal);
        alarm(m_seconds);
        do {
          m_user.op();
          m_samples++;
        } while (!m_gotSignal);
#endif
      }
    } else {
      m_seconds = 0;
    }

    m_user.sync();
    clock_t time_start = clock();
    for (i = 0; i < m_samples; i++) m_user.op();
    m_user.sync();
    clock_t time_end = clock();
    m_period = (float) (time_end - time_start) / (float) CLOCKS_PER_SEC;

    // std::cout << m_samples << " iterations, " << m_period <<" seconds"
    //           << std::endl;
  }

  /*!
   */
  void printResults()
  {
    fflush(stdout);
    if (!m_headerPrinted && m_printHeader) printHeader();
    m_headerPrinted = true;
    int count = (int) ((float) m_samples * m_factor);
    std::cout << std::setw(m_nameLength) << std::left
              << m_name.substr(0, m_nameLength).c_str() << " "
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

  Bench_user & getBenchUser() { return m_user; }
    
private:
  Bench_user m_user;

  std::string m_name;
  bool m_printHeader;
  bool m_headerPrinted;
  float m_factor;
  int m_seconds;
  int m_samples;
  float m_period;
  int m_iterations;
};

CGAL_END_NAMESPACE

#endif
