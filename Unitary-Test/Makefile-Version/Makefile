# ==============================================================
# Template Makefile for Eigen 3
# ==============================================================

CPP = g++
CPPFLAGS = -O2
INCLUDE_DIR = -I/usr/local/include/eigen3
LINK_LIB = 
PROGRAM_NAME = solver

$(PROGRAM_NAME): $(PROGRAM_NAME).cpp
	$(CPP) $(INCLUDE_DIR) $(CPPFLAGS) -o $@ $< $(LINK_LIB)

clean:
	rm -r $(PROGRAM_NAME)