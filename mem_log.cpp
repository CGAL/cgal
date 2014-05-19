#include <Windows.h>
#include <psapi.h>

#include <iostream>
#include <fstream>
//#include <ios>
#include <mem_log.h>
using namespace std;
void LogMyMemoryUsage()
{
    PPROCESS_MEMORY_COUNTERS pMemCountr = new PROCESS_MEMORY_COUNTERS;

    if(GetProcessMemoryInfo(GetCurrentProcess(),pMemCountr, 
            sizeof(PROCESS_MEMORY_COUNTERS)))
    {
       // TDateTime t;

        std::string csMsg;
         //"Time = " + t.CurrentTime().TimeString() + 
		//csMsg =    " – Memory Usage(kb) = " + std::string(((pMemCountr->PagefileUsage/1024))) + "\n";
		double usage = pMemCountr->PagefileUsage/1024;
		cout << "Mem usage : "<< usage << endl;
  /*      ofstream myfile;
        myfile.open ("C:\\memorylog.txt", ios::out | ios::app);
		myfile << csMsg.c_str();
        myfile.close();*/
    }
    //lastMemoryUsage = (pMemCountr->PagefileUsage/1024);
    delete pMemCountr;
}