# $Id$

############# MACROS ##############
.SECONDEXPANSION:

include ../config.mk

# ******************************
NAME    = optimizer
INCDIR  = \
	  -I../utils \
	  -I../math \
	  -I../graph \
	  -I../geometry \
	  -I../data
SOURCES = $(wildcard *.cc)

######### COMMON MAKEFILE #########
include ../rules.mk

############# TARGETS #############
all: depend $(LIBRARY)

# **************************
# create dependencies
# **************************
depend: $(DEPS)

# and include them
ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif

$(LIBRARY): $(OBJECTS)
	rm -f $@
	ar rv $@ $^
	ranlib $@

.PHONY: clean
clean:
	rm -f  *.$(OBJ) $(LIBRARY) $(DEPS)

.PHONY: tags
tags:
	ctags *.h *.cc
