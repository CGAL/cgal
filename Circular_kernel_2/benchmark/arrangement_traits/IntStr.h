
#include <string>
#include <sstream>

template <class T>
bool from_string(T &t,
                 const std::string &s,
                 std::ios_base & (*f)(std::ios_base&))
{
   std::istringstream iss(s);
   return !(iss>>f>>t).fail();
}


template <class T>
std::string to_string(T t, std::ios_base & (*f)(std::ios_base&)=std::dec)
{
   std::ostringstream oss;
   oss << f << t;
   return oss.str();
}


