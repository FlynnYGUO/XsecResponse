# -*- mode: makefile -*- 
#
# Define some directories and search paths
#
HERE  := $(PWD)
LIBDIR = $(HERE)/lib
OBJDIR = $(HERE)/.obj
DEPDIR = $(HERE)/.dep
EXEDIR = $(HERE)/bin

LOCAL_INC = -I$(HERE)/../src_lib -I$(HERE)/../src_bin -I/vols/dune/nk3717/duneanaobj_build_v2CAFs//build/install/include 
vpath % $(HERE)/../src_lib
vpath % $(HERE)/../src_bin
ROOTVERSION = $(shell root-config --version)
#
# Set up compiler and flags
#
CC        = g++
LD_FLAGS  = $(shell root-config --cflags) -fPIC -L$(LIBDIR) -Wl,-undefined,dynamic_lookup $(LOCAL_INC)
CC_FLAGS  = $(OPT) -Wall -Werror $(shell root-config --cflags) -fPIC $(LOCAL_INC)

ROOT_LIBS  = $(shell root-config --libs) -lMinuit

LIB_FLAGS = $(LD_FLAGS) $(ROOT_LIBS) -shared 
EXE_FLAGS = $(LD_FLAGS) $(ROOT_LIBS) $(LOCAL_LIBS) -Iinclude/

#
# Debug or Optimized
#
OPT = -g
#OPT = -O2


#
# Define library and executable source files
#
#LIB        = $(LIBDIR)/libXsecVary1D_NDGAr.so
LIB        = $(LIBDIR)/libXsecVary2D_NDGAr.so
SRCS       = $(wildcard ../src_lib/*.cc)
EXESRCS    = $(wildcard ../src_bin/*.cc)

# Make sure the executables can find the library above
#LOCAL_LIBS = -L$(LIBDIR) -lXsecVary1D_NDGAr -L/vols/dune/nk3717/duneanaobj_build/build/install/lib -lduneanaobj_StandardRecord
LOCAL_LIBS = -L$(LIBDIR) -lXsecVary2D_NDGAr -L/vols/dune/nk3717/duneanaobj_build_v2CAFs/build/install/lib -lduneanaobj_StandardRecord



#
# Automatically create library object and executable names
#
LIBOBJS = $(addprefix $(OBJDIR)/,$(notdir $(SRCS:.cc=.o)))
EXES = $(addprefix $(EXEDIR)/,$(notdir $(EXESRCS:.cc=)))
df = $(DEPDIR)/$(*F)


#
# Build Rules
#

all : directories $(LIB) $(EXES) 

$(EXEDIR)/% : $(OBJDIR)/%.o  $(LIB)
	@echo "Linking executable ..."
	$(CC) $(EXE_FLAGS) -o $@ $< $(STATIC_LIBS)

$(LIB) : $(LIBOBJS)
	@echo "Linking shared library ..."
	$(CC) --shared $(LIB_FLAGS) -o $@ $^

# Looks complicated -- used to automatically generate
# the dependency files as part of the build process.
# Probably overkill for this simple structure, but 
# it can't hurt.
$(OBJDIR)/%.o : %.cc
	@echo "Compiling object ..."
	$(CC) -c $(CC_FLAGS) -MD -MF $(df).d -o $@ $<
	@cp $(df).d $(df).P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $(df).d >> $(df).P; \
	    rm -f $(df).d

-include $(SRCS:%.cc=$(DEPDIR)/%.P)

#
# Not using dictionaries now, but may need them in the future
#
$(EXEDIR)/%Dict.C: %.cc Linkdef.h
	@echo "Generating dictionary ..."
	rootcint -f $@ -c $(CC_FLAGS) -p $^


directories : $(LIBDIR) $(OBJDIR) $(DEPDIR) $(EXEDIR)
$(LIBDIR) : 
	mkdir -p $@
$(OBJDIR) : 
	mkdir -p $@
$(DEPDIR) :
	mkdir -p $@
$(EXEDIR) :
	mkdir -p $@



.PHONY : clean
clean : 
	@rm -f $(DEPDIR)/*
	@rm -f $(LIBOBJS)
	@rm -f $(DICT)
	@rm -f $(DICT:.C=.h)
	@rm -f $(LIB)
	@rm -f $(EXES)
