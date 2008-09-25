// Author: Laurent Saboret

#ifndef RedirectIOToFile_H
#define RedirectIOToFile_H

#include <stdio.h>


// This function opens a file and points stdout and stderr to it.
// You must close it later with fclose().
FILE* RedirectIOToFile(const char* log_file_name)
{
    FILE* log_file;
    
    // Open file
    log_file = fopen(log_file_name, "w");
    if (log_file == NULL)
      return NULL;

    // redirect unbuffered STDOUT to the log file
    *stdout = *log_file;
    setvbuf( stdout, NULL, _IONBF, 0 );

    // redirect unbuffered STDERR to the log file
    *stderr = *log_file;
    setvbuf( stderr, NULL, _IONBF, 0 );	
    
    return log_file;
}


#endif // RedirectIOToFile_H
