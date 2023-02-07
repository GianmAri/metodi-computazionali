CXX		:= `root-config --cxx`
CXXFLAGS	:= -std=c++11 -O3 -Wall `root-config --cflags`
LDFLAGS		:= `root-config --ldflags`
SRC		:=./src
INC		:=./include
INCLUDES	:=-I./$(INC) -I`root-config --incdir`
OBJ		:=./obj/
LIB		:=-lm `root-config --glibs`

MAIN		:= $(OBJ)/main.o
CLASSES		:= particle
EADERS		:= $(CLASSES:%=$(INC)/%.h)
SOURCES		:= $(CLASSES:%=$(SRC)/%.cxx)
OBJECTS		:= $(CLASSES:%=$(OBJ)/%.o)
PROJECTNAME	:= project

default: main

$(OBJ)%.o: $(SRC)/%.cxx
	@echo Compiling $< ...
	@if ! [ -d $(OBJ) ] ; then mkdir -pv $(OBJ); fi
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

main: $(OBJECTS) $(MAIN)
	@echo Linking $^ to $@
	$(CXX)  $^ -o $@ $(LDFLAGS) $(LIB)

clean:
	@rm -rfv $(OBJ)
	@rm -fv main
	@rm -fv $(PROJECTNAME).zip

zippa:
	@zip -r $(PROJECTNAME).zip src include Makefile
