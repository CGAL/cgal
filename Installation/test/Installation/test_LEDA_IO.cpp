#include <LEDA/numbers/rational.h>

#include <cstdlib>
#include <iostream>
#include <sstream>

const char* tests[] = { "0\n", "1/2\n", "3 4", "0", "1/2" };

int main()
{
  for(int i = 0, end = sizeof(tests)/sizeof(char*);
      i < end; ++i)
  {
    std::cout << "input: \"" << tests[i] << "\"\n";
    std::stringstream input(tests[i]); // no \n: EOF instead
    leda::rational a;
    input >> a;
    if(!input.fail()) {
      std::cout << "a is " << a << std::endl;
      if(std::cout.fail()) {
        std::cerr << "std::count.fail() is set!\n";
        return EXIT_FAILURE;
      }
    } else {
      std::cerr << "input.fail() is set!\n";
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
