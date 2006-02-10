 
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

// int main()
// {
//    int i;
//    float f;
//    // the third parameter of from_string() should be 
//    // one of std::hex, std::dec or std::oct
//    if(from_string<int>(i, std::string("ff"), std::hex)){
//       std::cout<<i<<std::endl;
//    }
//    else{
//       std::cout<<"from_string failed"<<std::endl;
//    }
//    if(from_string<float>(f, 
//                                std::string("123.456"),
//                                std::dec))
//    {
//       std::cout<<f<<std::endl;
//    }
//    else{
//       std::cout<<"from_string failed"<<std::endl;
//    }
//    return 0;
// } 

/* output:
255
123.456
*/

template <class T>
std::string to_string(T t, std::ios_base & (*f)(std::ios_base&)=std::dec)
{
   std::ostringstream oss;
   oss << f << t;
   return oss.str();
}

// int main()
// {
//    // the second parameter of to_string() should be one of 
//    // std::hex, std::dec or std::oct
//    std::cout<<to_string<long>(123456, std::hex)<<std::endl;
//    std::cout<<to_string<long>(123456, std::oct)<<std::endl;
//    return 0;
// } 

 
