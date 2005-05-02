#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "taucs_structure.h"

#define INPUT_LEN 1024
char  reply[INPUT_LEN];

char* ostype;
char* variant      = 0;
int   build_flags;
char  pathsep      = '/';
int   win32        = 0;

void mark_dependants()
{
  int i,j,k;
  int go = 1;
  while (go) {
    go = 0;
    for (i=0; modules[i].name; i++) {
      if (modules[i].state == 1) {
	for (j=0; modules[i].depends[j]; j++) {
	  for (k=0; modules[k].name; k++) {
	    if (!strcmp(modules[i].depends[j],modules[k].name) && modules[k].state!=1) {
	      go=1;
	      modules[k].state=1;
	    }
	  }
	}
      }
    }
  }
}

void user_input()
{
  int i;

  printf("\n");

  for (i=0; modules[i].name; i++)
    modules[i].state = -1;

  for (i=0; modules[i].name; i++) {
    if (modules[i].state != 1) { /* could be 1 due to a dependency */

      printf("module %s? (default is %s) ",
	     modules[i].name,
	     modules[i].inc_exc ? "yes" : "no");
      fgets(reply, INPUT_LEN, stdin);
      
      while (modules[i].state == -1) {
	switch (reply[0]) {
	case 'y':
	case 'Y':
	  modules[i].state = 1;
	  break;
	case 'n':
	case 'N':
	  modules[i].state = 0;
	  break;
	case '\r':
	case '\n':
	  modules[i].state = modules[i].inc_exc;
	  break;
	default:
	  printf("<%c> is not a valid response, please type 'y' or 'n' ",reply[0]);
	  fgets(reply, INPUT_LEN, stdin);
	}
      }
    }
    mark_dependants();
  }
}

void map_sources_to_targets()
{
  int i,j,k,target;

  for (i=0; files[i].name; i++)
    files[i].state = 0;

  for (i=0; modules[i].name; i++) {
    if ( ! modules[i].state ) continue; /* if we don't make the module, skip              */
    if ( ! modules[i].target ) continue; /* if it has no target, don't bother with sources */
    /*printf(" including module %s\n",modules[i].name);*/
    for (k=0; files[k].name; k++)
      if (!strcmp(modules[i].target,files[k].name)
	  &&
	  (files[k].flags & target_types)
	  ) {
	target=k;
	/*printf("<<< %s/%s\n",modules[i].name,files[k].name);*/
	files[k].state = -1;
      }

    for (j=0; modules[i].sources[j]; j++)
      for (k=0; files[k].name; k++)
	if ((!strcmp(modules[i].sources[j],files[k].name)) 
	    &&
	    (files[k].flags & source_types)) {
	  /*printf(">>> %s/%s\n",modules[i].name,files[k].name);*/
	  files[k].state=target;
	}
  }
}

/*
char* new_get_path(char* path)
{
  char* p;

  p = (char*) malloc(strlen(path)+1);
  strcpy(p,path);

  for (; *p; p++)
    if (*p == '/' || *p == '\\')
      *p = (!strcmp(platform,"win32")) ? '\\' : '/';
    
  return path;
}
*/
  
char* get_extension(int flags)
{
  /*int i,j;*/

  if (flags & object    ) return "$(OBJEXT)";
  if (!(flags & object) && (flags & fcsource)) return "$(F2CEXT)";
  if (flags & hsource   ) return ".h";
  if (flags & csource   ) return ".c";
  if (flags & cilksource) return ".c";
  if (flags & lib       ) return "$(LIBEXT)";
  if (flags & executable) return "$(EXEEXT)";

  /*
  for (i=0; extensions[i].name; i++) {
    if (extensions[i].flags & flags) {
      for (j=0; extensions[i].platforms[j].name; j++)
	if (!strcmp(extensions[i].platforms[j].name,platform))
	  return extensions[i].platforms[j].extension;
    }
  }
  */
  fprintf(stderr,"error: internal error, missing extension for flag=%x\n",flags);
  exit(1);
  
  return 0;
}

char* get_full_name(int i,int modifier) {
  char* name;
  int j;

  name = (char*) malloc(256);

  if (modifier & cilksource) {
    sprintf(name, "$(%s)%s%s",
	    files[i].path,
	    files[i].name,
	    ".cilk");
    return name;
  }

  if (modifier & object) {
    for (j=0; number_types[j].flag; j++) {
      if (number_types[j].flag & modifier) {
	sprintf(name, "$(%s)%s%s%s",
		"DIROBJ",
		files[i].name,
		number_types[j].extension,
		get_extension(object));
	return name;
      }
    }
    sprintf(name, "$(%s)%s%s",
	    "DIROBJ",
	    files[i].name,
	    get_extension(object));
    return name;
  }
  
  if (modifier & executable) {
    for (j=0; number_types[j].flag; j++) {
      if (number_types[j].flag & modifier) {
	sprintf(name, "$(%s)%s%s%s",
		"DIROBJ",
		files[i].name,
		number_types[j].extension,
		get_extension(executable));
	return name;
      }
    }
    sprintf(name, "$(%s)%s%s",
	    "DIROBJ",
	    files[i].name,
	    get_extension(executable));
    return name;
  }
  
  sprintf(name, "$(%s)%s%s",
	  files[i].path,
	  files[i].name,
	  get_extension(files[i].flags));
  return name;
}

char* get_lib_spec(int i) {
  char* name;
  /*return get_full_name(i,0);*/

  if (win32) {
    return get_full_name(i,0);
  } else {
    if (strncmp(files[i].name,"lib",3)) {
      fprintf(stderr,"error: library name must begin with 'lib' (%s)\n",files[i].name);
      exit(1);
    }
    name = (char*) malloc(256);
    sprintf(name, "-L$(%s) -l%s",
	    files[i].path,
	    files[i].name + 3);
  }

  return name;
}

void emit_configfile(char* configuration_name)
{
  int i;
  FILE* f;
  char  name[256];

  /* create the configuration directory */

  if (variant)
    sprintf(name,"%s%c%s%s",configdir,pathsep,ostype,variant);
  else
    sprintf(name,"%s%c%s",  configdir,pathsep,ostype);

  if (win32) {
    mkdir(configdir);
    mkdir(name);
  }
  else {
    mkdir(configdir,0777);
    mkdir(name,0777);
  }

  /* now create the configuration file */

  if (variant)
    sprintf(name,"%s%c%s%s%c%s",configdir,pathsep,ostype,variant,pathsep,configfile);
  else
    sprintf(name,"%s%c%s%c%s",  configdir,pathsep,ostype,pathsep,configfile);
  
  f = fopen(name,"w");
  if (!f) {
    fprintf(stderr,"error: could not open config file <%s>\n",name);
    exit(1);
  }

  /*
  f = fopen(configfile,"w");
  if (!f) {
    fprintf(stderr,"error: could not open config file <%s>\n",configfile);
    exit(1);
  }
  */
  fprintf(f, "/* This is an automatically generated file */\n");
  fprintf(f, "/* Configuration name: %s */\n",configuration_name);
  fprintf(f,"#define %s_OSTYPE %s\n",prefix,ostype);
  fprintf(f,"#define %s_VARIANT %s\n",prefix,variant?variant:"none");
  fprintf(f,"#define OSTYPE_%s\n",ostype);
  fprintf(f,"#define OSTYPE_VARIANT_%s\n",variant?variant:"none");
  for (i=0; modules[i].name; i++) {
    if ( ! modules[i].state  ) continue;
    fprintf(f,"#define %s_CONFIG_%s\n",prefix,modules[i].name);
  }
  fprintf(f, "\n");
  
  fclose(f);
}  


void emit_makefile(char* configuration_name)
{
  int i,j,k;
  FILE* f;
  char* base;
  char  name[256];

  /* create the configuration directory */

  if (variant)
    sprintf(name,"%s%c%s%s",configdir,pathsep,ostype,variant);
  else
    sprintf(name,"%s%c%s",  configdir,pathsep,ostype);

  if (win32) {
    mkdir(configdir);
    mkdir(name);
  }
  else {
    mkdir(configdir,0777);
    mkdir(name,0777);
  }

  /* now create the makefile */

  if (variant)
    sprintf(name,"%s%c%s%s%c%s",configdir,pathsep,ostype,variant,pathsep,makefile);
  else
    sprintf(name,"%s%c%s%c%s",  configdir,pathsep,ostype,pathsep,makefile);

  /* create an "alias" */
  f = fopen("makefile","w");
  if (!f) {
    fprintf(stderr,"error: could not open makefile <%s>\n",makefile);
    exit(1);
  }
  fprintf(f,"include %s\n",name);
  fclose(f);

  /* now create the real makefile */
  f = fopen(name,"w");
  if (!f) {
    fprintf(stderr,"error: could not open makefile <%s>\n",name);
    exit(1);
  }

  fprintf(f, "OSTYPE = %s\n", ostype );  
  if (variant)
    fprintf(f, "OSTYPE_VARIANT = %s\n", variant );  

  fprintf(f, "CONFIGURATION = %s\n", configuration_name);

  /*
  fprintf(f, "\n");
  if (!strcmp(platform,"win32")) {
    fprintf(f, "EXEOUT=/Fe\n" );
    fprintf(f, "OBJOUT=/Fe\n" );
    fprintf(f, "DEFFLG=/D\n"  );
    fprintf(f, "EXEEXT=.exe\n"  );
  } else {
    fprintf(f, "EXEOUT=-o\n" );
    fprintf(f, "OBJOUT=-o\n" );
    fprintf(f, "DEFFLG=-D\n"  );
    fprintf(f, "EXEEXT=\n"  );
  }
  */

  /*
  fprintf(f, "\n");
  for (i=0; paths[i].name; i++) {
    fprintf(f,"%s=%s\n",paths[i].name,get_path(paths[i].name));
  }
  */

  fprintf(f, "INCS = ");
  /*fprintf(f, "\\\n  %s",configfile);*/
  fprintf(f, "\\\n  %s%c%s%s%c%s",
	  configdir,pathsep,ostype,
	  variant?variant:"",pathsep,configfile);
  fprintf(f, "\\\n  %s%c%s%s%c%s",
	  configdir,pathsep,ostype,
	  variant?variant:"",pathsep,testfile);
  for (i=0; files[i].name; i++) {
    base = files[i].name;
    if (files[i].flags & hsource)
      fprintf(f,"\\\n  %s",get_full_name(i,0));
  }
  fprintf(f, "\n");

  fprintf(f, "CILKC=$(CC)\n");
  fprintf(f, "CILKFLAGS=$(CFLAGS)\n");
  fprintf(f, "CILKOUTFLG=$(COUTFLG)\n");
  fprintf(f, "\n");

  fprintf(f, "\n");
  if (win32) {
    fprintf(f, "include config%cwin32.mk\n", pathsep);
    if (variant)
      fprintf(f, "include config%cwin32%s.mk\n", pathsep,variant);
  } else {
    fprintf(f, "include config%c$(OSTYPE).mk\n", pathsep);
    if (variant)
      fprintf(f, "include config%c$(OSTYPE)$(OSTYPE_VARIANT).mk\n", pathsep);
  }
  fprintf(f,"\n");

  fprintf(f, "LIBS =");
  for (i=0; modules[i].name; i++)
    if (modules[i].state) 
      for (k=0; modules[i].extralib[k]; k++)
	fprintf(f," $(%s)",modules[i].extralib[k]);
  fprintf(f," $(%s)",f77_lib);
  fprintf(f," $(%s)",c_lib);
  fprintf(f,"\n");

  fprintf(f, "default: all\n" );  
  fprintf(f, "\n");

  fprintf(f, "include config%cstd.mk\n", pathsep);
  fprintf(f, "\n");

  fprintf(f, "%s%c%s%s%c%s: $(DIROBJ)exists.log",
	  configdir,pathsep,ostype,
	  variant?variant:"",pathsep,testfile);
  for (i=0; files[i].name; i++) {
    if (files[i].flags & test)
      fprintf(f, " %s",get_full_name(i,0));
  }
  fprintf(f, "\n");

  for (i=0; files[i].name; i++) {
    if ( ! (files[i].flags & test) ) {
      continue;
    }
    
    if (files[i].flags & csource) {
      fprintf(f, "\t- $(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \\\n");
      fprintf(f, "\t  %s \\\n",get_full_name(i,0));
      fprintf(f, "\t  $(COUTFLG)%s\n", get_full_name(i,object));
    }

    if (files[i].flags & cilksource) {
      fprintf(f, "\t- $(CILKC) -c $(CILKFLAGS) $(STDDEFS) $(STDINCS) \\\n");
      fprintf(f, "\t  %s \\\n",get_full_name(i,0));
      fprintf(f, "\t  $(CILKOUTFLG)%s\n", get_full_name(i,object));
    }
	 
    fprintf(f,"\t- $(LD) $(LDFLAGS) \\\n");
    fprintf(f,"\t  $(LOUTFLG)%s \\\n",get_full_name(i,executable));
    fprintf(f,"\t  %s $(LIBS)\n",get_full_name(i,object));
    
    fprintf(f,"\t- %s",get_full_name(i,executable));
    fprintf(f, " %s%c%s%s%c%s\n",
	    configdir,pathsep,ostype,
	    variant?variant:"",pathsep,testfile);
  }

  for (i=0; files[i].name; i++) {
    base = files[i].name;

    if ( ! files[i].state ) {
      /*printf(">>> state=0 for %s\n",base);*/
      continue;
    }

    if ( (files[i].flags & csource) 
	 || (files[i].flags & cilksource) ) {
      int generated = 0;
      for (j=0; number_types[j].flag; j++) {
	if ( ! (number_types[j].flag & build_flags) ) continue;
	if (files[i].flags & number_types[j].flag) {
	  generated = 1;

	  /*if ((build_flags & cilksource) && (files[i].flags & cilksource)) {*/
	  if (files[i].flags & cilksource) {
	    fprintf(f, "%s: %s $(INCS) $(STDDEPS)\n",
		    get_full_name(i,object | number_types[j].flag),
		    get_full_name(i,0));
	    /*	    
	    fprintf(f, "\t/bin/cp %s %s\n",get_full_name(i,0),get_full_name(i,cilksource));
	    fprintf(f, "\t$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \\\n");
	    fprintf(f, "\t-D%s \\\n",number_types[j].define);
	    fprintf(f, "\t-DTAUCS_CORE_CILK \\\n");
	    fprintf(f, "\t%s \\\n",get_full_name(i,cilksource));
	    fprintf(f, "\t$(COUTFLG)%s\n",get_full_name(i,object | number_types[j].flag));
	    fprintf(f, "\t/bin/rm %s\n",get_full_name(i,cilksource));
	    */
	    fprintf(f, "\t$(CILKC) -c $(CILKFLAGS) $(STDDEFS) $(STDINCS) \\\n");
	    fprintf(f, "\t-D%s \\\n",number_types[j].define);
	    fprintf(f, "\t%s \\\n",get_full_name(i,0));
	    fprintf(f, "\t$(CILKOUTFLG)%s\n",get_full_name(i,object | number_types[j].flag));
	  } else {
	    fprintf(f, "%s: %s $(INCS) $(STDDEPS)\n",
		    get_full_name(i,object | number_types[j].flag),
		    get_full_name(i,0));
	    fprintf(f, "\t$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \\\n");
	    fprintf(f, "\t-D%s \\\n",number_types[j].define);
	    fprintf(f, "\t%s \\\n",get_full_name(i,0));
	    fprintf(f, "\t$(COUTFLG)%s\n",get_full_name(i,object | number_types[j].flag));
	  }

	}
      }
      if (!generated) {
	/*if ((build_flags & cilksource) && (files[i].flags & cilksource)) {*/
	if (files[i].flags & cilksource) {
	  fprintf(f, "%s: %s $(INCS) $(STDDEPS)\n",
		  get_full_name(i,object),
		  get_full_name(i,0));
	  /*
	  fprintf(f, "\t/bin/cp %s %s\n",get_full_name(i,0),get_full_name(i,cilksource));
	  fprintf(f, "\t$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \\\n");
	  fprintf(f, "\t-DTAUCS_CORE_CILK \\\n");
	  fprintf(f, "\t%s \\\n",get_full_name(i,cilksource));
	  fprintf(f, "\t$(COUTFLG)%s\n",
		  get_full_name(i,object));
	  fprintf(f, "\t/bin/rm %s\n",get_full_name(i,cilksource));
	  */
	  fprintf(f, "\t$(CILKC) -c $(CILKFLAGS) $(STDDEFS) $(STDINCS) \\\n");
	  fprintf(f, "\t%s \\\n",get_full_name(i,0));
	  fprintf(f, "\t$(CILKOUTFLG)%s\n",
		  get_full_name(i,object));
	} else {
	  fprintf(f, "%s: %s $(INCS) $(STDDEPS)\n",
		  get_full_name(i,object),
		  get_full_name(i,0));
	  fprintf(f, "\t$(CC) -c $(CFLAGS) $(STDDEFS) $(STDINCS) \\\n");
	  fprintf(f, "\t %s \\\n",get_full_name(i,0));
	  fprintf(f, "\t$(COUTFLG)%s\n",
		  get_full_name(i,object));
	}

      }
    }

    if (files[i].flags & fcsource) {
      fprintf(f, "%s: %s $(STDDEPS)\n",
	      get_full_name(i,object),
	      get_full_name(i,0));
      fprintf(f, "\t$(FC) -c $(FFLAGS) \\\n");
      fprintf(f, "\t%s \\\n",get_full_name(i,0));
      fprintf(f, "\t$(FOUTFLG)%s\n",get_full_name(i,object));
    }
  }

  for (i=0; files[i].name; i++) {
    base = files[i].name;

    if ( ! files[i].state ) {
      /*printf(">>> state=0 for %s\n",base);*/
      continue;
    }

    if ( (files[i].flags & lib) ) {
      fprintf(f, "%s_content = ",base);
      for (k=0; files[k].name; k++) {
	int   generated = 0;
	if (files[k].state != i) continue; /* this source does not contribute to this lib */
	for (j=0; number_types[j].flag; j++) {
	  if ( ! (number_types[j].flag & build_flags) ) continue;
	  if (files[k].flags & number_types[j].flag) {
	    generated = 1;
	    fprintf(f, "\\\n  %s ",get_full_name(k,object | number_types[j].flag));
	  }
	}
	if (!generated) {
	  fprintf(f, "\\\n  %s ",get_full_name(k,object));
	}
      }
      fprintf(f,"\n");

      fprintf(f, "%s: $(%s_content) $(STDDEPS)\n",get_full_name(i,0),base);
      fprintf(f,"\t- $(RM) %s\n",get_full_name(i,0));
      fprintf(f,"\t$(AR) $(AOUTFLG)%s $(%s_content)\n",get_full_name(i,0),base);
      fprintf(f,"\t$(RANLIB) %s\n",get_full_name(i,0));
    }

    if ( (files[i].flags & executable) ) {
      /*
      fprintf(f, "%s_content = ",base);
      for (k=0; files[k].name; k++) {
	if (files[k].state != i) continue;
	fprintf(f, "\\\n  %s ",get_full_name(k,object));
      }
      fprintf(f,"\n");
      */

      /*fprintf(f, "%s: $(%s_deps) $(%s)\n",get_full_name(i,0),base,files[i].path);*/

      fprintf(f, "%s: $(STDDEPS)",get_full_name(i,0));

      for (k=0; files[k].name; k++) {
	if (files[k].state != i) continue;
	fprintf(f, " %s",get_full_name(k,object));
      }
      for (k=0; files[k].name; k++) {
	if ( ! (files[k].flags & lib) ) continue; /* we assume an exe depends on all libs  */
	if ( ! (files[k].state) )       continue; /* except if the lib is not built at all */
	fprintf(f, " %s",get_full_name(k,0));
      }
      fprintf(f,"\n");

      fprintf(f,"\t$(LD) $(LDFLAGS) \\\n");
      fprintf(f,"\t$(LOUTFLG)%s \\\n",get_full_name(i,0));
      for (k=0; files[k].name; k++) {
	if (files[k].state != i) continue;
	fprintf(f, "\t%s \\\n",get_full_name(k,object));
      }
      /*fprintf(f,"\t$(%s_content) \\\n",base);*/
      for (k=0; files[k].name; k++) {
	if ( ! (files[k].flags & lib) ) continue; /* we assume an exe depends on all libs  */
	if ( ! (files[k].state) )       continue; /* except if the lib is not built at all */
	fprintf(f, "\t%s \\\n",get_lib_spec(k));
      }
      fprintf(f,"\t$(LIBS)\n");
    }
  }

  /* now emit the list of targets */

  fprintf(f,"all:");
  fprintf(f," $(STD_PRE_TARGETS)");
  for (i=0; files[i].name; i++) {
    if ( ! files[i].state ) continue;
    if (files[i].flags & target_types) {
      fprintf(f," %s",get_full_name(i,0));
    }
  }
  fprintf(f,"\n");

  fclose(f);
}

void parse_modulename(char* s)
{
  int   i;
  int   ok = 0;
  int   include_exclude;
  char  module_name[64];
  char* p;
  char* q = module_name;

  if (!s || !(*s)) {
    fprintf(stderr,"error: module name is missing\n");
    exit(1);
  }

  if (*s == '!') {
    include_exclude=0;
    p = s+1;
  } else {
    include_exclude=1;
    p = s;
  }

  while (*p && ((q-module_name)<64) && *p != ' ' && *p != '\t' && *p != '\r' && *p != '\n') {
    *q = *p;
    p++;
    q++;
  }
  *q = 0;
  
  for (i=0; modules[i].name; i++) {
    if (!strcmp(modules[i].name,module_name)) {
      modules[i].state=include_exclude;
      ok=1;
      break;
    }
  }
  
  if (!ok) {
    fprintf(stderr,"error: module name <%s> is not a recognized module\n",module_name);
    exit(1);
  }
}

void parse_infile(char* name)
{
  FILE* f;
  int   i;
  int   next;
  char  s[256];

  f = fopen(name,"r");
  if (!f) {
    fprintf(stderr,"error: could not open input file <%s>\n",name);
    exit(1);
  }

  for (i=0; modules[i].name; i++)
    modules[i].state = -1;

  sprintf(s,"%s_CONFIG_BEGIN",prefix);
  next=0;
  while(!feof(f)) {
    fgets(reply, INPUT_LEN, f);
    if (!strncmp(reply,s,strlen(s))) { next=1; break; }
  }
  if (!next) {
    fprintf(stderr,"error: could not open find label <%s> in input file <%s>\n",s,name);
    exit(1);
  }

  next=0;
  while(!feof(f)) {
    fgets(reply, INPUT_LEN, f);

    sprintf(s,"%s_CONFIG_END",prefix);
    if (!strncmp(reply,s,strlen(s))) { next=1; break; }

    sprintf(s,"%s_CONFIG_DEFAULT ",prefix);
    if (!strncmp(reply,s,strlen(s))) {
      int ok = 0;
      if (!strncmp(reply+strlen(s),"OFF",3)) {
	for (i=0; modules[i].name; i++) modules[i].state=0;
	ok=1;
      } 
      if (!strncmp(reply+strlen(s),"ON",2)) {
	for (i=0; modules[i].name; i++) modules[i].state=1;
	ok=1;
      }
      if (!ok) {
	fprintf(stderr,"error: label <%s> in input file <%s> not followed by ON or OFF\n",
		s,name);
	exit(1);
      }
    }

    sprintf(s,"%s_CONFIG ",prefix);
    if (!strncmp(reply,s,strlen(s))) parse_modulename(reply+strlen(s));

  }
  if (!next) {
    sprintf(s,"%s_CONFIG_END",prefix);
    fprintf(stderr,"error: could not open find label <%s> in input file <%s>\n",s,name);
    exit(1);
  }

  for (i=0; modules[i].name; i++)
    sprintf(s,"%s_CONFIG %s\n",prefix, modules[i].name);
  
  sprintf(s,"%s_CONFIG_END\n",prefix);

  fclose(f);
}

void emit_outfile(char* name)
{
  FILE* f;
  int i;

  f = fopen(name,"w");
  if (!f) {
    fprintf(stderr,"error: could not open output file <%s>\n",name);
    exit(1);
  }

  fprintf(f,"%s_CONFIG_BEGIN\n",prefix);
  fprintf(f,"%s_CONFIG_DEFAULT OFF\n",prefix);

  for (i=0; modules[i].name; i++)
    if (modules[i].state) fprintf(f,"%s_CONFIG %s\n",prefix, modules[i].name);
    else                  fprintf(f,"%s_CONFIG !%s\n",prefix, modules[i].name);
  
  fprintf(f,"%s_CONFIG_END\n",prefix);

  fclose(f);
}

int arg_get_boolean(char* arg, char* pattern, int* x) {
  if (!strcmp(arg,pattern)) {
    *x = 1;
    return 1;
  }

  if (arg[0] && arg[0]=='!' && !strcmp(arg+1,pattern)) {
    *x = 0;
    return 1;
  }

  return 0;
}

int arg_get_string(char* arg, char* pattern, char** x) {
  if (!strncmp(arg,pattern,strlen(pattern))) {
    *x = arg + strlen(pattern);
    return 1;
  }
  
  return 0;
}

char* get_configuration_name(char* fname)
{
  int found_prefix;
  char* name;

  if (!fname) return "anonymous";

  name = (char*) malloc(strlen(fname)+1);
  
  found_prefix=0;
  if (strrchr(fname,'/')) {
    strcpy(name,strrchr(fname,'/')+1);
    found_prefix=1;
  }
  if (strrchr(fname,'\\')) {
    strcpy(name,strrchr(fname,'\\')+1);
    found_prefix=1;
  }
  if (!found_prefix) strcpy(name,fname);

  if (strchr(name,'.')) *( strchr(name,'.') ) = 0;

  return name;
}

int main(int argc, char* argv[])
{
  int i,j;

  int   interactive  = 0;
  char* outfile      = 0;
  char* defaultfile  = 0;
  char* module       = 0;

  /* we set up the defaults */

  for (j=0; modules[j].name; j++)
    modules[j].state = modules[j].inc_exc;

  for (i=1; i<argc; i++) {

    /*fprintf(stderr,"argv[%d]\n",i);*/

    if (arg_get_boolean(argv[i],"help"       ,&j /*dummy*/ )) {
      fprintf(stderr,"  \n");
      fprintf(stderr,"usage: %s [options]\n",argv[0]);
      fprintf(stderr,"\n");
      fprintf(stderr,"  This program builds a makefile and a configuration\n");
      fprintf(stderr,"  include file for %s according to a selection\n",prefix);
      fprintf(stderr,"  of modules.\n");
      fprintf(stderr,"  With no options, it builds a configuration according\n");
      fprintf(stderr,"  to built-in default rules. The options are:\n");
      fprintf(stderr,"  variant=<name>        an OSTYPE variant\n");
      /*
      fprintf(stderr,"  makefile=<filename>   overrides the built-in name (%s)\n",makefile);
      fprintf(stderr,"  configfile=<filename> overrides the built-in name (%s)\n",configfile);
      */
      fprintf(stderr,"  win32                 builds a makefile for nmake\n");
      fprintf(stderr,"                        (otherwise to unix/GNU makes)\n");
      fprintf(stderr,"  interactive           asks the users which modules to select\n");
      fprintf(stderr,"  in=<filename>         module selection taken from a file\n");
      fprintf(stderr,"  out=<filename>        module selection stored to a file\n");
      fprintf(stderr,"  module=<modulename>   selects a module\n");
      fprintf(stderr,"  module=!<modulename>  de-selects a module\n");
      fprintf(stderr,"  \n");
      fprintf(stderr,"  There can be several module= options, and they should\n");
      fprintf(stderr,"  normally follow in= options; the module= options can\n");
      fprintf(stderr,"  then override selections that the file makes.\n");
      fprintf(stderr,"  \n");
      fprintf(stderr,"  Configuration files have the following structure:\n");
      fprintf(stderr,"  ...\n");
      fprintf(stderr,"  %s_CONFIG_BEGIN\n",prefix);
      fprintf(stderr,"  %s_CONFIG_DEFAULT [ON/OFF]\n",prefix);
      fprintf(stderr,"  %s_CONFIG_MODULE SOME_MODULE\n",prefix);
      fprintf(stderr,"  %s_CONFIG_MODULE !ANOTHER_MODULE\n",prefix);
      fprintf(stderr,"  %s_CONFIG_MODULE A_THIRD_MODULE\n",prefix);
      fprintf(stderr,"  %s_CONFIG_END\n",prefix);
      fprintf(stderr,"  ...\n");
      fprintf(stderr,"  These lines can be embedded in C comments, etc, but\n");
      fprintf(stderr,"  they must begin at the beginnig of lines. The second\n");
      fprintf(stderr,"  line in the example allows you to select or deselct\n");
      fprintf(stderr,"  all the modules together.\n");
      fprintf(stderr,"  \n");
      fprintf(stderr,"  The list of modules (and whether they are included in\n");
      fprintf(stderr,"  the default configuration is:\n");
      fprintf(stderr,"  \n");
      for (j=0; modules[j].name; j++)
      fprintf(stderr,"  %-20s %s\n",modules[j].name,
	     modules[j].inc_exc ? "included" : "excluded");
      fprintf(stderr,"  \n");
      exit(0);
     }

    if (arg_get_string (argv[i],"variant="   ,&variant     )) continue;
    /*
    if (arg_get_string (argv[i],"makefile="  ,&makefile    )) continue;
    if (arg_get_string (argv[i],"configfile=",&configfile  )) continue;
    */
    if (arg_get_string (argv[i],"out="       ,&outfile     )) continue;

    if (arg_get_boolean(argv[i],"win32"      ,&win32       )) {
      if (win32) pathsep = '\\';
      continue;
    }

    if (arg_get_boolean(argv[i],"interactive",&interactive )) {
      continue;
    }

    if (arg_get_string (argv[i],"in="        ,&defaultfile )) {
      parse_infile(defaultfile);
      mark_dependants();
      continue;
    }

    if (arg_get_string (argv[i],"module="    ,&module      )) {
      parse_modulename(module);
      mark_dependants();
      continue;
    }

    fprintf(stderr,"error: command line argument <%s> is invalid\n",argv[i]);
    exit(1);
  }

  ostype = getenv("OSTYPE");
  if (!ostype) {
    fprintf(stderr,"error: the environment variable OSTYPE is not set\n");
    exit(1);
  }

  if (interactive) user_input();

  for (i=0; modules[i].name; i++)
    if (modules[i].state)
      build_flags |= modules[i].build_flags;

  map_sources_to_targets();

  emit_configfile(get_configuration_name(defaultfile));

  emit_makefile(get_configuration_name(defaultfile));

  if (outfile)
    emit_outfile(outfile);

  printf("%s%s\n",ostype,variant?variant:"");

  return 0;
}
