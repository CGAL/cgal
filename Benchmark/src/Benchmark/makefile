# Copyright (c) 1997  Tel-Aviv University (Israel).
# All rights reserved.
#
# This file is part of CGAL (www.cgal.org); you may redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; version 2.1 of the License.
# See the file LICENSE.LGPL distributed with CGAL.
#
# This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
# WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# $URL$
# $Id$
# 
# Author(s)     : Efi Fogel <efif@post.tau.ac.il>

CXX ?=g++
CXXFLAGS ?=-Wall
CPPFLAGS = -I. -I../../include -g
AR ?=ar
LEX =flex
BISON ?=bison
INSTALL ?=install
RANLIB ?=ranlib
PREFIX ?=/usr/local

TARGET =libCGALBenchmark.a
CPPSOURCES =Benchmark_visitor.cpp
CPPSOURCES+=Bench.cpp
CPPSOURCES+=Bench_option_parser.cpp
BISONSOURCES =benchmark_parser.y
LEXSOURCES =benchmark_lexer.l
DIRT =benchmark_parser.cpp benchmark_parser.hpp benchmark_parser.output
DIRT+=benchmark_lexer.cpp

CPPOBJECTS =$(CPPSOURCES:.cpp=.o)
BISONOBJECTS =$(BISONSOURCES:.y=.o)
LEXOBJECTS =$(LEXSOURCES:.l=.o)
OBJECTS = $(CPPOBJECTS) $(BISONOBJECTS) $(LEXOBJECTS)

default: all
all: $(TARGET)
$(TARGET): $(OBJECTS)
	@echo Creating library $@ ...
	$(AR) crus $@ $?

# install: ${TARGET}
# 	${INSTALL} -c -d ${PREFIX}/lib
# 	${INSTALL} -c -p -m 644 ${TARGET} ${PREFIX}/lib/${TARGET}
# 	${RANLIB} ${PREFIX}/lib/${TARGET}

clean:
	-rm -f *.o $(TARGET) $(DIRT)

.SUFFIXES:
.SUFFIXES: .o .l .y .cpp

%.cpp %.hpp : %.y
	$(BISON) -d $< -o $*.cpp

%.cpp : %.l
	$(LEX) -8 -o$@ $<

Benchmark_visitor.o: benchmark_parser.hpp
benchmark_lexer.o: benchmark_parser.hpp
