PROG =	Hausdorff.out

SRCS =  

OBJS =  constants.o initialise.o hausdorff.o deriv.o blas.f lbfgsb.o linpack.o timer.o optimization.o rotations.o simulate.o main.o

### gfortran options
LIBS =
#FC = gfortran
#FFLAGS = -g -ffree-form -Wall -Wextra -O2 -fbounds-check 
LDFLAGS= 
RUNCPP = 
   
### ifort options   
FC = ifort
FFLAGS = -g -O -free -ftz -ip -ipo -qopenmp -parallel -prec-div -prec-sqrt 

all: $(PROG)

$(PROG): $(OBJS) $(EXTRAS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(EXTRAS) $(LIBS) 


clean:
	rm $(PROG)
	rm *.o

