# makefile 
FC = mpiifort
CHECK = -CB -g -traceback -check all,noarg_temp_created -debug all
FCFLAGS = -r8 -O3 -free -mcmodel=large -heap-arrays 10 -shared-intel -fp-model precise -qmkl
ifdef DEBUG
  FCFLAGS += $(CHECK)
endif
FINCLUDE = -I${HOME}/.local/include -I/opt/libs-intel-oneapi/netcdf-4.7.4/include
LDLIBS = -L${HOME}/.local/lib -lfftw3 -L/opt/libs-intel-oneapi/netcdf-4.7.4/lib -lnetcdff


# code paths
VPATH = ./src/

# objects

LIST = cal_ens.F 
a.out: cal_ens.o 
cal_ens.o: cal_ens.F

LIST_o = $(LIST:.F=.o)
target = a.out 


all: $(target)

$(LIST_o): %.o: %.F
	$(FC) $(FCFLAGS) $(FINCLUDE) -c $<

$(target) : $(LIST_o)
	$(FC) $(FCFLAGS) $(FINCLUDE) $^ -o $@ $(LDLIBS)

clean:
	rm -rf *.o *.mod ace_model.exe


