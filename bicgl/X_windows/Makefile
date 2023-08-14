def: all

include ../../BIC_PL/Makefile.include

LIBS = -lX11 $(BIC_PL_LIBS)
LINT_LIBS =  $(BIC_PL_LINTLIBS)

OPT = $(OPT_g)

OPT_g = -g
OPT_O = -O

LIB_NAME = x_windows
XWINDOW_LIB      = lib$(LIB_NAME).a
XWINDOW_LINT_LIB = llib-l$(LIB_NAME).ln

OBJECTS = x_windows.o x_draw.o

INCLUDE = $(BIC_PL_INCLUDE) -IInclude -I../Include

#-----------------------------------------------------------------------------

X_WINDOW_PROTOTYPE_FILE = Include/x_window_prototypes.h

$(X_WINDOW_PROTOTYPE_FILE): $(OBJECTS:.o=.c)
	@echo "#ifndef  DEF_X_WINDOW_PROTOTYPES"         >  $@
	@echo "#define  DEF_X_WINDOW_PROTOTYPES"         >> $@
	@extract_functions - $(OBJECTS:.o=.c) >> $@
	@echo "#endif"                                      >> $@

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

$(XWINDOW_LIB): $(X_WINDOW_PROTOTYPE_FILE) $(OBJECTS)
	@if( -e $@ ) \rm -f $@
	$(MAKE_LIBRARY) $@ $(OBJECTS)
	$(RANLIB) $@

$(XWINDOW_LINT_LIB): $(OBJECTS:.o=.ln)
	@echo "--- Linting ---"
	$(LINT) -x -u -o $(LIB_NAME) $(OBJECTS:.o=.ln)

all: $(XWINDOW_LIB) $(XWINDOW_LINT_LIB) test_x_windows

test_x_windows: $(XWINDOW_LIB) $(XWINDOW_LINT_LIB) \
                 test_x_windows.o test_x_windows.ln
	$(CC) test_x_windows.o -o $@ $(XWINDOW_LIB) $(LIBS)
	$(LINT) -x -u test_x_windows.ln $(XWINDOW_LINTLIB) $(LINT_LIBS)
