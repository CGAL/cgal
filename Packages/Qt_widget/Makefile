# This makefile needs a file named ./.env
# This .env file can be used to declare environment variables.  It is
# usefull to have different environment per working copy of the
# repository. It use that two work with two different version of Qt on
# two working copies.
# For example, I have two different .env file:
# one for Qt-2.3:
#>export QTDIR=/usr/lib/qt2
#>export CGAL_MAKEFILE=/opt/CGAL/make/makefile_i686_Linux-2.4.8-26mdk_g++-2.96.
# and one for Qt-3.0:
#>export QTDIR=/opt/local/qt-3
#>export CGAL_MAKEFILE=/opt/CGAL/make/makefile_i686_Linux-2.4.8-26mdk_g++-3.0.1-3.0.1
# If that bothers you, just comment the following line:
# [Laurent Rineau 2002/01/30]

include .env

DEMO_DIRS= \
./demo/Qt_widget/Partition_2 \
./demo/Qt_widget \
./demo/Qt_widget/Triangulation_2 \
./demo/Qt_widget/Convex_Hull_2 \
./demo/Qt_widget/basic/tutorial1 \
./demo/Qt_widget/basic/tutorial2 \
./demo/Qt_widget/basic/tutorial3 \
./demo/Qt_widget/basic/tutorial4 \
./demo/Qt_widget/basic/tutorial5 \
./demo/Qt_widget/basic/tutorial6 \
./demo/Qt_widget/basic/tutorial7 \
./demo/Qt_widget/basic/tutorial8

# make the lib and all demos
all:
	make -C src/CGALQt -f makefile_Qt
	for dir in ${DEMO_DIRS}; do make -C $$dir; done

# clean all object files
clean: 
	make -C src/CGALQt -f makefile_Qt clean
	for dir in ${DEMO_DIRS}; do make -C $$dir clean; done

# clean moc files too
clean-moc: clean
	find -name '*.moc' -exec rm \{\} \;

# make update: update the working directory
update:
	@cvs -q update -d -P

# make diff: diff between the cvs repository and the working copy
diff:
	@cvs -q diff | less
