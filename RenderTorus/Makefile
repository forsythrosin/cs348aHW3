include makefile.in

INCLUDE = -I$(OPENMESH_INCLUDE_DIR) -Iinclude/
CPPFLAGS = -O3 -fPIC 
LDFLAGS = -O3 -lGL -lGLU
LIB = -lglut -lOpenMeshCored -lOpenMeshToolsd -Wl,-rpath,$(OPENMESH_LIB_DIR)
TARGET = renderTorus
OBJS = objs/main-torus.o

default: $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -L$(OPENMESH_LIB_DIR) $(LIB) -o $(TARGET)
	
objs/main-torus.o: src/main-torus.cpp
	$(CPP) -c $(CPPFLAGS) src/main-torus.cpp -o objs/main-torus.o $(INCLUDE)

clean:
	rm $(OBJS) $(TARGET) -f
