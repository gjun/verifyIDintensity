# version information
TOOL = verifyIDintensity
# 0.0.1 - Initial version
VERSION = 0.0.1
SOURCES = verifyIDintensity.cpp
#OPTFLAG?=-ggdb -O0 -ggdb
OPTFLAG?=-O2 

# default installation directory
#INSTALLDIR=/usr/local/bin

# default C++ compiler
CXX=g++ 

DATE=$(shell date)
NODE=$(shell uname -n)
#
# default compilation flags are 
#
# CFLAGS=-02 -I./libsrc/ -I./$(TOOL)
# 
# The following special options may also be added to the default
# 
#      Option                        Effect
#      -D_FILE_OFFSET_BITS=64        Enables support for swapfiles larger
#             	                      than 2GB on supported systems (tested
#                                    on Linux)
#      -D__USE_LONG_INT              Enables support for markers with up
#                                    to 64 alleles (default is 32). Tested
#                                    on systems where gcc supports the long
#                                    long data type and on Windows.
# 

INCLUDES=-I. -I../include

CFLAGS=-pipe -Wall -Wno-trigraphs -fno-rtti $(OPTFLAG) $(INCLUDES) -D__ZLIB_AVAILABLE__ -D__STDC_LIMIT_MACROS 

# executable file names and locations
TARGET = $(TOOL)

# Source File Set
# For best results, consider editing this manually ...
TOOLBASE = 
TOOLHDR = $(TOOLBASE:=.h)
TOOLONLYSRC = $(TOOLBASE:=.cpp)
TOOLSRC = $(TOOLBASE:=.cpp) $(SOURCES)
TOOLOBJ = $(TOOLSRC:.cpp=.o)
TOOLONLYOBJ = $(TOOLONLYSRC:.cpp=.o)


# Utility Library File Set
LIBRARY=
NEWLIBRARY=

# private parameters
#FETCHDIR=$(HOME)/code
#DISTRIBDIR=$(HOME)/code/distrib/$(TOOL)-$(VERSION)

# make everything
all : $(TARGET) $(NEWLIBRARY)

# dependencies for executables
$(TARGET) : $(LIBRARY) $(TOOLOBJ) $(BINDIR)
	$(CXX) $(CFLAGS) -o $@ $(TOOLOBJ) $(LIBRARY) -lm -lz

$(TOOLOBJ): $(TOOLHDR) $(LIBHDR)


$(NEWLIBRARY) : $(TOOLONLYOBJ)
	ar -cr $@ $(TOOLONLYOBJ)
	ranlib $@

clean :
	-rm -f *.o $(TARGET) *~ $(NEWLIBRARY)

.c.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.c 

.cpp.o : 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DVERSION="\"$(VERSION)\""



.SUFFIXES : .cpp .c .o .X.o $(SUFFIXES)

DFLAGS=-Y

cleandepend:
	        makedepend -- $(DFLAGS) --

depend:
	        makedepend -- $(DFLAGS) -- $(TOOLSRC) >/dev/null 2>&1

# DO NOT DELETE THIS LINE -- make depend depends on it

verifyIDintensity.o: common.h MarkerStat.h Stat.h Intensity.h SampleStat.h
