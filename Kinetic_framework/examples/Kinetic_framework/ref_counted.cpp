#include <iostream>
#include <CGAL/Kinetic/Ref_counted.h>
struct Foo: CGAL::Kinetic::Ref_counted<Foo>
{

  ~Foo() {
    std::cout << "Foo: bye, bye." << std::endl;
  }
};

struct Bar: CGAL::Kinetic::Ref_counted<Bar>
{

  Bar(Foo::Handle fp): fp_(fp){}

  ~Bar() {
    std::cout << "Bar: bye, bye." << std::endl;
  }
  Foo::Handle fp_;
};

int main(int, char*[])
{

  Foo::Handle fp= new Foo;
  std::cout << "Clearing fp.\n";
  fp= NULL;
  fp = new Foo;
  Bar::Handle bp= new Bar(fp);
  std::cout << "Clearing fp again.\n";
  fp=NULL;
  std::cout << "Clearing bp.\n";
  bp=NULL;

  return EXIT_SUCCESS;
}
