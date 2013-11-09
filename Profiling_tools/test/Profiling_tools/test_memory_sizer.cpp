// Test program for CGAL::Memory_sizer.
//
// Sylvain Pion

#include <CGAL/Memory_sizer.h>
#include <iostream>

int main()
{
    typedef CGAL::Memory_sizer::size_type size_type;

    CGAL::Memory_sizer mem;

    size_type v0 = mem.virtual_size();
    size_type r0 = mem.resident_size();

    std::cout << "Virtual  size = " << v0
              << " (" << (v0>>20) << " MB)" << std::endl;
    std::cout << "Resident size = " << r0
              << " (" << (r0>>20) << " MB)" << std::endl;

    std::cout << "Now allocating 64MB" << std::endl;

    char *p = new char[64*1024*1024];

    size_type v1 = mem.virtual_size();
    size_type r1 = mem.resident_size();

    std::cout << "Virtual  size = " << v1
              << " (" << (v1>>20) << " MB)" << std::endl;
    std::cout << "Resident size = " << r1
              << " (" << (r1>>20) << " MB)" << std::endl;

    std::cout << "Now touching it" << std::endl;

    for (int i=0; i<64*1024; ++i)
        p[1024*i] = 1;

    size_type v2 = mem.virtual_size();
    size_type r2 = mem.resident_size();

    std::cout << "Virtual  size = " << v2
              << " (" << (v2>>20) << " MB)" << std::endl;
    std::cout << "Resident size = " << r2
              << " (" << (r2>>20) << " MB)" << std::endl;

    delete [] p;
    return 0;
}
