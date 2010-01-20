

#include <boost/thread.hpp>

class MessageQueue
{
public:
  typedef boost::function1<void,int> Fnc;
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
    
//    if ( status_ == PREUNACTIVATED )
//      myCondition.notify_all();
    
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
pooledThread(MessageQueue* queue, int i)
{
  MessageQueue::Receive_type f = queue->receive(true);
  while ( f.first )
  {
    f.second(i);
    f = queue->receive(false);
  }
  
  std::stringstream tmp;
  tmp << "end of thread " << i << "\n";
  std::cout << tmp.str();
} 