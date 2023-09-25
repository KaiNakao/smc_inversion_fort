COMPILER_F90 = mpifrtpx
COMPILER_CPP = mpiFCCpx
# FLAGS_F90   = -Kfast,openmp -SSL2
# FLAGS_CPP   = -Kfast,openmp -SSL2
FLAGS_F90   = -g -Kopenmp -SSL2 -Haesofux
FLAGS_CPP   = -g -Kopenmp

# CFLAGS   = -Wall -Wextra -Wno-sign-compare -g -fopenmp -fsanitize=address,undefined 
# FFLAGS   = -g -fopenmp
TARGET   = ./main
OBJDIR   = ./obj
SOURCES_F90 = $(wildcard src/*.f90)
OBJECTS_F90 = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES_F90:.f90=.o)))
SOURCES_CPP = $(wildcard src/*.cpp)
OBJECTS_CPP = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES_CPP:.cpp=.o)))
OBJECT_DC3D = $(addprefix $(OBJDIR)/, DC3Dfortran.o) 

$(TARGET): $(OBJECTS_F90) $(OBJECT_DC3D) $(OBJECTS_CPP)
	$(COMPILER_F90) -o $@ $^ $(FLAGS_F90) --linkstl=libfjc++

obj/sort.o: src/sort.cpp
	$(COMPILER_CPP) -o $@ -c $^ $(FLAGS_CPP) 

obj/init.o: src/init.f90
	$(COMPILER_F90) -o $@ -c $^ $(FLAGS_F90) 

obj/gfunc.o: src/gfunc.f90 obj/DC3Dfortran.o
	$(COMPILER_F90) -o $@ -c $^ $(FLAGS_F90) 

obj/smc_slip.o: src/smc_slip.f90
	$(COMPILER_F90) -o $@ -c $^ $(FLAGS_F90) 

obj/smc_fault.o: src/smc_fault.f90 obj/init.o obj/gfunc.o obj/smc_slip.o obj/sort.o
	$(COMPILER_F90) -o $@ -c $^ $(FLAGS_F90) 
	
obj/main.o: src/main.f90 obj/init.o obj/gfunc.o obj/smc_fault.o
	$(COMPILER_F90) -o $@ -c $^ $(FLAGS_F90) 

obj/DC3Dfortran.o: src/DC3Dfortran.f
	$(COMPILER_F90) -o $@ -c $^ $(FLAGS_F90) 

all: clean $(TARGET)

clean:
	rm -f obj/* $(TARGET) *.mod src/*.mod *.f90
