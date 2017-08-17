INCLUDE_DIR = include -I/usr/local/include/eigen3
SOURCE_DIR = src
OBJECT_DIR = obj
CC = g++
DEFFLAGS = -DVTK
CFLAGS = -Wall -Wno-write-strings -O2 -I$(INCLUDE_DIR)
CC_SOURCES = $(wildcard src/*)
STRING_OBJ_AUX = $(CC_SOURCES:.cpp=.o)
STRING_OBJ = $(subst src/,,$(STRING_OBJ_AUX))
CC_OBJ = $(patsubst %,$(OBJECT_DIR)/%,$(STRING_OBJ))
LINKER_LIBS = -lm 
PROGRAM_NAME = heat2D

all: $(PROGRAM_NAME)

$(PROGRAM_NAME): $(CC_OBJ)
	$(CC) $(DEFFLAGS) $(CFLAGS) -o $@ $^ $(LINKER_LIBS)

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	$(CC) $(DEFFLAGS) $(CFLAGS) -o $@ -c $< $(LINKER_LIBS)

clean:
	rm -f $(OBJECT_DIR)/*.o $(PROGRAM_NAME)
	rm -f VTK/*

plot:
	cd Data; python plot.py

rerun:
	rm -f VTK/*.vtk

print-%  : ; @echo $* = $($*)
