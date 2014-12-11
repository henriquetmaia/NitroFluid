#######################################################################################################

# Specify library locations here (add or remove "#" marks to comment/uncomment lines for your platform)

# Mac OS X
DDG_INCLUDE_PATH      = -I/usr/local/include/ -I/opt/local/include/
DDG_LIBRARY_PATH      = -L/usr/local/lib/ -L/opt/local/lib/
DDG_BLAS_LIBS         = -framework Accelerate
DDG_SUITESPARSE_LIBS  = -lspqr -lumfpack -lcholmod -lmetis -lcolamd -lccolamd -lcamd -lamd -lm -lsuitesparseconfig
DDG_OPENGL_LIBS       = -framework OpenGL -framework GLUT

# # Linux
# DDG_INCLUDE_PATH      =
# DDG_LIBRARY_PATH      =
# DDG_BLAS_LIBS         = -llapack -lblas -lgfortran
# DDG_SUITESPARSE_LIBS  = -lspqr -lcholmod -lmetis -lcolamd -lccolamd -lcamd -lamd -lm
# DDG_OPENGL_LIBS       = -lglut -lGL -lGLU -lX11

# # Windows / Cygwin
# DDG_INCLUDE_PATH      = -I/usr/include/opengl -I/usr/include/suitesparse
# DDG_LIBRARY_PATH      = -L/usr/lib/w32api -L/usr/lib/suitesparse
# DDG_BLAS_LIBS         = -llapack -lblas
# DDG_SUITESPARSE_LIBS  = -lspqr -lcholmod -lcolamd -lccolamd -lcamd -lamd -lm
# DDG_OPENGL_LIBS       = -lglut32 -lglu32 -lopengl32

#######################################################################################################

TARGET = ddg
CC = g++
LD = g++
CFLAGS = -O3 -Wall -Wno-deprecated -Werror -pedantic  $(DDG_INCLUDE_PATH) -I./include -I./src -DNDEBUG
LFLAGS = -O3 -Wall -Wno-deprecated -Werror -pedantic $(DDG_LIBRARY_PATH) -DNDEBUG
LIBS = $(DDG_OPENGL_LIBS) $(DDG_SUITESPARSE_LIBS) $(DDG_BLAS_LIBS)

## !! Do not edit below this line -- dependencies can be updated by running ./update ##################

OBJS = obj/Camera.o obj/Complex.o obj/Image.o obj/LinearEquation.o obj/Quaternion.o obj/Real.o obj/Shader.o obj/Variable.o obj/Vector.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(LD) $(LFLAGS) $(OBJS) $(LIBS) -o $(TARGET)

obj/Camera.o: src/Camera.cpp include/Camera.h include/Quaternion.h include/Vector.h include/Complex.h 
	$(CC) $(CFLAGS) -c src/Camera.cpp -o obj/Camera.o

obj/Complex.o: src/Complex.cpp include/Complex.h 
	$(CC) $(CFLAGS) -c src/Complex.cpp -o obj/Complex.o

obj/Image.o: src/Image.cpp include/Image.h 
	$(CC) $(CFLAGS) -c src/Image.cpp -o obj/Image.o

obj/LinearEquation.o: src/LinearEquation.cpp include/LinearEquation.h include/LinearPolynomial.h include/Variable.h 
	$(CC) $(CFLAGS) -c src/LinearEquation.cpp -o obj/LinearEquation.o

obj/Quaternion.o: src/Quaternion.cpp include/Quaternion.h include/Vector.h include/Complex.h 
	$(CC) $(CFLAGS) -c src/Quaternion.cpp -o obj/Quaternion.o

obj/Real.o: src/Real.cpp include/Real.h 
	$(CC) $(CFLAGS) -c src/Real.cpp -o obj/Real.o

obj/Shader.o: src/Shader.cpp include/Shader.h 
	$(CC) $(CFLAGS) -c src/Shader.cpp -o obj/Shader.o

obj/Variable.o: src/Variable.cpp include/Variable.h 
	$(CC) $(CFLAGS) -c src/Variable.cpp -o obj/Variable.o

obj/Vector.o: src/Vector.cpp include/Vector.h 
	$(CC) $(CFLAGS) -c src/Vector.cpp -o obj/Vector.o


clean:
	rm -f $(OBJS)
	rm -f $(TARGET)
	rm -f $(TARGET).exe

