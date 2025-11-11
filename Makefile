MODULE = nucphys
SRC1 = utils/woodsax.f90
SRC2 = utils/fnucutils.f90
SRC3 = utils/simf.f90
OUTDIR = utils/
TARGET = $(OUTDIR)$(MODULE).so
# 
.PHONY: all clean
# 
all: $(TARGET)
# 
$(TARGET): $(SRC)
	f2py -m $(MODULE) -c $(SRC1) $(SRC2) $(SRC3)
	mv $(MODULE)*.so $(OUTDIR)/

 
clean:
	rm -f $(OUTDIR)$(MODULE).cpython-*.so 
