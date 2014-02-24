include makefile.in

INCLUDE = -I$(OPENMESH_INCLUDE_DIR) -Iinclude/ -I$(EIGEN_DIR)
CPPFLAGS = -O3 -fPIC -DEIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS -DEIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET 
LDFLAGS = -O3 -lGL -lGLU
LIB = -lglut -lOpenMeshCored -lOpenMeshToolsd -Wl,-rpath,$(OPENMESH_LIB_DIR)
TARGET = drawMesh
OBJS = objs/main-curvature.o objs/curvature.o

default: $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -L$(OPENMESH_LIB_DIR) $(LIB) -o $(TARGET)
	
objs/main-curvature.o: src/main-curvature.cpp
	$(CPP) -c $(CPPFLAGS) src/main-curvature.cpp -o objs/main-curvature.o $(INCLUDE)

objs/curvature.o: src/curvature.cpp
	$(CPP) -c $(CPPFLAGS) src/curvature.cpp -o objs/curvature.o $(INCLUDE)

clean:
	rm $(OBJS) $(TARGET) -f
