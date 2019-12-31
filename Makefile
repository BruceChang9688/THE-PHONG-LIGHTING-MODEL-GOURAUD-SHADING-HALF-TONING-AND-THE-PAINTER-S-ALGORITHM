CC = /usr/bin/g++

#OPENGL LIBS FOR LINUX
#GLLIB :=  -lGL -lGLEW -lGLU -lglut
#OPENGL LIBS FOR MAC
#$(CC) $(CFLAGS) $^  -o $@ /usr/local/lib/libglfw.3.3.dylib /usr/local/lib/libGLEW.2.1.0.dylib
GLLIB := -framework OpenGL -framework GLUT

#COMPILER FLAGS
CCFLAGS :=

#include directories
#should include gl.h glut.h etc...
INCDIR := -I/usr/local/include
LDLIBS := $(GLLIB)

TARGET = Project3
OBJS = main.o


all: $(TARGET)


$(TARGET): $(OBJS)
	$(CC)  $^ $(CCFLAGS) $(LDLIBS)  -o $@

%.o : %.cpp
	$(CC) $(CCFLAGS) -o $@ -c $(INCDIR) $<

clean:
	rm -rf $(OBJS) $(TARGET)

