MODULE = fwoodsax
SRC = utils/woodsax.f90
OUTDIR = utils/
TARGET = $(OUTDIR)$(MODULE).so

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	f2py -m $(MODULE) -c $(SRC)
	mv $(MODULE)*.so $(OUTDIR)/

clean:
	rm -f $(OUTDIR)$(MODULE).cpython-*.so 
