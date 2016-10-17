#include <stdio.h>
#include <CGAL/exceptions.h>
#include <boost/thread.hpp>

std::string format(std::string,bool);

enum {
  PRINT_ERROR_MSG=-1
};

class MessageQueue
{
public:
  typedef boost::function2<void,unsigned int,std::string> Fnc;
  typedef std::pair<bool,Fnc> Receive_type;
  
  MessageQueue(): status_(ACTIVATED), running_nb_(0) {}
  
  void send( Fnc const& f)
  {
    boost::mutex::scoped_lock lock( myMutex );
    myQueue.push_back( f );
    myCondition.notify_one();
  }
  
  Receive_type receive(bool first)
  {
    boost::mutex::scoped_lock lock( myMutex );
    if (!first) --running_nb_;
    
    while ( myQueue.empty() )
    {
      if (status_ == PREUNACTIVATED && running_nb_ == 0)
      { 
        status_ = UNACTIVATED;
        myCondition.notify_all();
      }
      
      if (status_ == UNACTIVATED)
        return std::make_pair(false,Fnc());

      myCondition.wait( lock );
    }
    
    std::pair<bool,Fnc> result (true,myQueue.front());
    myQueue.pop_front();
    ++running_nb_;
    return result;
  }
  
  void deactivate()
  { 
    boost::mutex::scoped_lock lock( myMutex );
    status_ = PREUNACTIVATED;
    myCondition.notify_one();
  }
  
private:
  enum Activation { ACTIVATED, UNACTIVATED, PREUNACTIVATED };
  
  boost::mutex myMutex;
  boost::condition_variable myCondition;
  std::deque< Fnc > myQueue;
  Activation status_;
  int running_nb_;
  int threads_nb_;
};



void
pooledThread(MessageQueue* queue, unsigned int i)
{
  MessageQueue::Receive_type f = queue->receive(true);
  while ( f.first )
  {
    try
    {
      f.second(i,std::string());
    }
    catch (CGAL::Failure_exception e)
    {
      std::cerr << "ERROR[" << i << "]:\n" << format(e.what(),false);
      f.second(PRINT_ERROR_MSG, std::string("\nERROR:\n"));
      f.second(PRINT_ERROR_MSG, format(e.what(),false));
    
#ifdef MESHER_TESTER_USE_CGAL_EXCEPTION_WITH_STACK    
      std::cerr << format(e.stack(),true) << std::endl;
      f.second(PRINT_ERROR_MSG, format(e.stack(),true));
#endif
    }
    catch (std::logic_error e)
    {
      std::cerr << "ERROR[" << i << "]:\n" << format(e.what(),false);
      f.second(PRINT_ERROR_MSG, std::string("\nERROR:\n"));
      f.second(PRINT_ERROR_MSG, format(e.what(),false));
    }
    catch (...)
    {
      std::cerr << "ERROR[" << i << "]: Unknown error\n";
      f.second(PRINT_ERROR_MSG, std::string("ERROR: Unknown error\n"));
    }

    std::stringstream msg;
    msg << "[" << i << "] End of job\n";
    std::cout << msg.str();
    f = queue->receive(false);
  }
  
  std::stringstream tmp;
  tmp << "end of thread " << i << "\n";
  std::cout << tmp.str();
}


std::string format(std::string input, bool remove_template)
{
  std::stringstream result;
 
  // Remove <>
  std::size_t pos = input.find('\074');
  while ( remove_template && pos != std::string::npos )
  {
    size_t init = pos;
    int nb_open=1;
    int nb_close=0;

    ++pos;
    while( pos < input.length() && nb_close != nb_open )
    {
      switch ( input.at(pos) )
      {
      case '\074': ++nb_open; break;
      case '\076': ++nb_close; break;
      }
      ++pos;
    }
    input.erase(init,pos-init);
    pos = input.find('\074',init+1);
  }

  // Add '\n' at the end of the string
  if ( input.at(input.size()-1) != '\n' )
    input.push_back('\n');

  // Add "  ! " at the begining of each line
  size_t prev = 0;
  pos = input.find("\n");
  while ( pos != std::string::npos )
  {
    result << "  ! " << input.substr(prev, pos-prev) << "\n";
    prev = pos+1;
    pos = input.find("\n",prev);
  }

  return result.str();
}
