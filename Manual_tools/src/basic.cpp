/**************************************************************************
 
  basic.cpp
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : Assertions, indented output, etc.
  System    : C++ (g++)
  Author    : (c) 1995 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#include <basic.h>
#include <stdio.h>

extern "C" {
#include <unistd.h>
}


// Own assertion macro
// ================================================

void cc_assertion_error(const char *cond, const char *fname, int line) {
    cerr << endl << "fatal error: assertion '" << cond << "' failed in line " 
	 << line << " of file '" << fname << "'." << endl;
    abort();
}

// Indented output, using global indentation counter
// ====================================================

class Output_indentation {
    int n;
    int step;
public:
    Output_indentation( int stepSize = 4) : n(0), step( stepSize) {}
    void operator++(){ n += step;}
    void operator--(){ n -= step;}
    int indentation() { return n; }
};

ostream& operator<< (ostream& out, Output_indentation& ind) {
    const char* s = "                                "; // 32 spaces !!!
    int n = ind.indentation();
    while ( n > 32) { // 32 !!!
	n -= 32;
	out << s;
    }
    if ( n > 0)
	out << (s + 32 - n);
    return out;
}

static Output_indentation output_indentation;

ostream& indent( ostream& out){
    ++ output_indentation;
    return out;
}

ostream& outdent( ostream& out){
    -- output_indentation;
    return out;
}

ostream& ind_newline( ostream& out){
    out << endl << output_indentation;
    return out;
}

ostream& ind_space( ostream& out){
    out << output_indentation;
    return out;
}

int indentation_number() {
    return output_indentation.indentation();
}

// Substitute old style malloc, realloc, strdup ...
// ================================================
char* renew( char* old, size_t old_size, size_t new_size) {
    CC_Assert( old_size == 0 || old != 0);
    CC_Assert( new_size > old_size);
    char* cpy = new char[ new_size];
    if ( old && old_size > 0) {
	size_t min = ( old_size < new_size ? old_size : new_size);
	memcpy( cpy, old, min);
	delete[] old;
    }
    return cpy;
}

char* newstr( const char* src) {
    CC_Assert( src);
    if ( ! src)
        return 0;
    char* s = new char[ strlen( src) + 1];
    strcpy( s, src);
    return s;
}

int execute_shell_command( const string& cmd, std::ostream& out, std::ostream& err ) {  
  int  stdin_pipe[2];
  int  stdout_pipe[2];
  int  stderr_pipe[2];
  int  retval_pipe[2];
  
  char buffer[BUFSIZ+1];
  char buffer_err[BUFSIZ+1];
  
  int fork_result;
  int data_processed;
  int data_processed_err;
  
  if( pipe(stdin_pipe) == 0 &&
      pipe(stdout_pipe) == 0 &&
      pipe(stderr_pipe) == 0 &&
      pipe(retval_pipe) == 0 )
  {
          fork_result = fork();
          if(fork_result == -1) {
                  std::cerr << "Fork Failure" << std::endl;
                  exit(EXIT_FAILURE);
          } else if(fork_result == 0) {
                  /* Close the Child process' STDIN */
                  close(0);
  
                  /* Duplicate the Child's STDIN to the stdin_pipe file descriptor */
                  dup(stdin_pipe[0]);
  
                  /* Close the read and write to for the pipe for the child.  The child will now only be able to read from it's STDIN (which is our pipe). */ 
                  close(stdin_pipe[0]);
                  close(stdin_pipe[1]);
  
                  /* Close the Child process' STDOUT */
                  close(1);
                  dup(stdout_pipe[1]);
                  close(stdout_pipe[0]);
                  close(stdout_pipe[1]);
  
                  /* Close the Child process' STDERR */
                  close(2);
                  dup(stderr_pipe[1]);
                  close(stderr_pipe[0]);
                  close(stderr_pipe[1]);
                  int retval = system( cmd.c_str() );
                  write(retval_pipe[1], &retval, sizeof(retval) );
                  exit(EXIT_SUCCESS);
          } else { // parent process
                  /* Close STDIN for read & write and close STDERR for write */
                  close(stdin_pipe[0]);
                  close(stdin_pipe[1]);
                  close(stderr_pipe[1]);
                  while(1) {
                          data_processed_err=read(stderr_pipe[0],buffer_err,BUFSIZ);
                          err.write( buffer_err,data_processed_err );
                          if( data_processed_err == 0 ) break;
                  }
                  /* Close the read end of STDERR */
                  close(stderr_pipe[0]);
                  /* Close the write end of STDOUT */
                  close(stdout_pipe[1]);
  
                  while(1) {
                          data_processed=read(stdout_pipe[0],buffer,BUFSIZ);          
                          out.write( buffer, data_processed );
                          if(data_processed == 0) break;
                  }
                  close(stdout_pipe[0]);
                  int retval;
                  read(retval_pipe[0], &retval, sizeof(retval) );
                  return retval;
          }
  } else {
    std::cerr << "pipe error" << std::endl;
    exit(EXIT_FAILURE);
  }
  return -1;
}

// EOF //
