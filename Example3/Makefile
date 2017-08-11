# =====================================================
# Template Makefile for building Eigen 3 with Qt
# =====================================================

CPP = g++
CPPFLAGS = -Wall -O2
INCLUDE_DIR = -I/usr/local/include/eigen3 -I/usr/share/qt4/mkspecs/linux-g++-64 -I. -I/usr/include/qt4/QtCore -I/usr/include/qt4/QtGui -I/usr/include/qt4
LINK_LIB = -L/usr/lib/x86_64-linux-gnu -lQtGui -lQtCore -lpthread
PROGRAM_NAME = example

$(PROGRAM_NAME): $(PROGRAM_NAME).cpp
	$(CPP) $(INCLUDE_DIR) $(CPPFLAGS) -o $@ $< $(LINK_LIB)

clean:
	rm -r $(PROGRAM_NAME)