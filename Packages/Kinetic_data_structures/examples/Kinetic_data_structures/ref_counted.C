#include <iostream>
#include <CGAL/KDS/Ref_counted.h>
struct Foo: CGAL::KDS::Ref_counted<Foo> {

  ~Foo(){
    std::cout << "Foo: bye, bye." << std::endl;
  }
};


struct Bar: CGAL::KDS::Ref_counted<Bar> {
  
  Bar(Foo::Pointer fp): fp_(fp){}
  
  ~Bar(){
    std::cout << "Bar: bye, bye." << std::endl;
  }
  Foo::Pointer fp_;
};

int main(int, char*[]){

  Foo::Pointer fp= new Foo;
  std::cout << "Clearing fp.\n";
  fp= NULL;
  fp = new Foo;
  Bar::Pointer bp= new Bar(fp);
  std::cout << "Clearing fp again.\n";
  fp=NULL;
  std::cout << "Clearing bp.\n";
  bp=NULL;

  return EXIT_SUCCESS;
}
