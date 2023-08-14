def: all

include ../../BIC_PL/Makefile.include

OPT = $(OPT_g)

OPT_g = -g
OPT_O = -O

LIB_NAME = glx_windows
GLXWINDOW_LIB      = lib$(LIB_NAME).a
GLXWINDOW_LINT_LIB = llib-l$(LIB_NAME).ln

OBJECTS = glx_windows.o stored_font.o

INCLUDE = $(BIC_PL_INCLUDE) -IInclude -I../X_windows/Include -I../Include

#-----------------------------------------------------------------------------

GLX_WINDOW_PROTOTYPE_FILE = Include/glx_window_prototypes.h

$(GLX_WINDOW_PROTOTYPE_FILE): $(OBJECTS:.o=.c)
	@echo "#ifndef  DEF_GLX_WINDOW_PROTOTYPES"         >  $@
	@echo "#define  DEF_GLX_WINDOW_PROTOTYPES"         >> $@
	@extract_functions - $(OBJECTS:.o=.c) >> $@
	@echo "#endif"                                      >> $@

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

$(GLXWINDOW_LIB): $(GLX_WINDOW_PROTOTYPE_FILE) $(OBJECTS)
	@if( -e $@ ) \rm -f $@
	$(MAKE_LIBRARY) $@ $(OBJECTS)
	$(RANLIB) $@

$(GLXWINDOW_LINT_LIB): $(OBJECTS:.o=.ln)
	@echo "--- Linting ---"
	$(LINT) -x -u -o $(LIB_NAME) $(OBJECTS:.o=.ln)

all: $(GLXWINDOW_LIB) $(GLXWINDOW_LINT_LIB)
