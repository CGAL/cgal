#!/bin/sh
#
# TAUCS: this makes the test script for windows and unix
#

/bin/rm -f testscript.bat
echo -e "REM This is an automatically-generated script\r" > testscript.bat
echo -e "del testscript.log\r" >> testscript.bat
echo -e "echo TAUCS TEST LOG >  testscript.log\r" >> testscript.bat
echo -e "echo ============== >> testscript.log\r" >> testscript.bat
echo -e "echo Win32          >> testscript.log\r" >> testscript.bat
echo -e "echo ============== >> testscript.log\r" >> testscript.bat

/bin/rm -f testscript
echo -e "#!/bin/sh" > testscript
echo -e "### This is an automatically-generated script" >> testscript
echo -e "/bin/rm testscript.log" >> testscript
echo -e "echo 'TAUCS TEST LOG' > testscript.log" >> testscript
echo -e "hostname >> testscript.log" >> testscript
echo -e "uname    >> testscript.log" >> testscript
echo -e "date     >> testscript.log" >> testscript
echo -e "echo '==============' >> testscript.log" >> testscript
echo -e "echo 'trying to maximize stack size:' >> testscript.log" >> testscript
echo -e "ulimit -s >> testscript.log" >> testscript
echo -e "ulimit -s unlimited >> testscript.log" >> testscript
echo -e "ulimit -s >> testscript.log" >> testscript
echo -e "echo '==============' >> testscript.log" >> testscript

chmod 755 testscript

for f in progs/test_*.c ; do
  echo $f
  
  name=`basename $f .c`
  bs='\\'
  fs='/'
  echo -e "echo =============== >> testscript.log\r" >> testscript.bat
  echo -e "echo = $name = >> testscript.log\r"       >> testscript.bat
  echo -e "call configure in=progs${bs}$name.c %*\r" >> testscript.bat
#  echo -e "nmake /F build${bs}%TAUCS_LASTCONF%${bs}makefile clean\r" >> testscript.bat
  echo -e "nmake /F build${bs}%TAUCS_LASTCONF%${bs}makefile \r"      >> testscript.bat
  echo -e "bin${bs}%TAUCS_LASTCONF%${bs}$name >> testscript.log\r" >> testscript.bat
  echo -e "if errorlevel 1 goto :error_$name\r"   >> testscript.bat
  echo -e "echo = TEST PASSED ($name) >> testscript.log\r" >> testscript.bat
  echo -e "goto :next_$name\r"                    >> testscript.bat
  echo -e ":error_$name\r"                        >> testscript.bat
  echo -e "echo = TEST FAILED ($name) >> testscript.log\r" >> testscript.bat
  echo -e ":next_$name\r"                         >> testscript.bat
  echo -e "echo =============== >> testscript.log\r" >> testscript.bat

  echo -e "echo =============== >> testscript.log" >> testscript
  echo -e "echo = $name = >> testscript.log"       >> testscript
  echo -e ". ./configure in=progs${fs}$name.c \$*" >> testscript
  echo -e "echo last conf is \$TAUCS_LASTCONF >> testscript.log" >> testscript
#  echo -e "make -f build${fs}\${TAUCS_LASTCONF}${fs}makefile clean"  >> testscript
  echo -e "make -f build${fs}\${TAUCS_LASTCONF}${fs}makefile"        >> testscript
  echo -e "if bin${fs}\${TAUCS_LASTCONF}${fs}$name >> testscript.log ; then" >> testscript
  echo -e "echo = TEST PASSED $name >> testscript.log" >> testscript
  echo -e "else"                    >> testscript
  echo -e "echo = TEST FAILED $name >> testscript.log" >> testscript
  echo -e "fi"                         >> testscript
  echo -e "echo =============== >> testscript.log" >> testscript
done
