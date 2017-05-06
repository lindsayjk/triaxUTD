CC = mpicc
FC = mpif90
AR = ar
TARGET = triaxUTD.a

CFLAGS = -std=c++11
FFLAGS =
INC = -I$(TACC_GSL_INC) -I$(TACC_GSL_INC)/inc
C_SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp,%.o,$(C_SRCS)) mcmc_f90.o

all: $(TARGET)
$(TARGET): $(OBJS)
	$(AR) rc $@ $^
%.o: %.cpp
	$(CC) $(CFLAGS) $(INC) -c $<
mcmc_f90.o: mcmc.f90
	$(FC) $(FFLAGS) -c $< -o mcmc_f90.o
clean:
	rm -rf $(TARGET) *.o *.mod

.PHONY: all clean
