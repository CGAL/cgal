RM = del
.SILENT:

.SUFFIXES: .clean

pwd.exe : pwd.clean pwd.cpp
	$(CXX) pwd.cpp $(extralibs)
	move pwd.exe ..\bin 
	$(RM) pwd.obj
	if exist pwd.tds $(RM) pwd.tds

replinc.exe: replinc.clean replinc.c
	$(CXX) replinc.c 
	move replinc.exe ..\bin
	$(RM) replinc.obj 
	if exist replinc.tds $(RM) replinc.tds

all: pwd.exe

# all: pwd.exe replinc.exe

.IGNORE:


.cpp.clean:
	if exist $*.obj $(RM) $*.obj
	if exist $*.exe $(RM) $*.exe
	if exist $*.tds $(RM) $*.tds
	if exist ..\bin\$*.exe $(RM) ..\bin\$*.exe

.c.clean:
	if exist $*.obj $(RM) $*.obj
	if exist $*.exe $(RM) $*.exe
	if exist $*.tds $(RM) $*.tds
	if exist ..\bin\$*.exe $(RM) ..\bin\$*.exe

