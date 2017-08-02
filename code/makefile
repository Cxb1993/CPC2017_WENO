f77=mpif90 
opt= -cpp -DTEST -DSWAP -DCHECK -c -OPT:ieee_arith=1

srs=  parameters.f90  main.f90 \
Weno.f90       \
Interfaces.f90
OBJS=$(srs:.f90=.o)

%.o:%.f90
	$(f77) $(opt) -c $<

default: $(OBJS)
	$(f77) -O3  -o weno7.out $(OBJS)

clean:
	rm -f *.out *.o work.* *.pc *.mod *.rmod

