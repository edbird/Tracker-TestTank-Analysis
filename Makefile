

#CC = g++
#INC = -I. -I./Analysis 
#CFLAGS = -Wall -Wextra --std=c++11 -O3 `root-config --cflags` `root-config --libs` $(INC)
#LDFLAGS = `root-config --cflags` `root-config --libs` -pthread

#Compiler and Linker
CC := g++

#The Target Binary Program
TARGET := run.out

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR := src
INCDIR := include
BUILDDIR := obj
TARGETDIR := bin
#RESDIR := res
SRCEXT := cpp
DEPEXT := d
OBJEXT := o

#Flags, Libraries and Includes
#CFLAGS := -fopenmp -Wall -O3 -g
CFLAGS := -Wall -Wextra --std=c++11 -O3 `root-config --cflags` `root-config --libs` $(INC)
#LIB := -fopenmp -lm -larmadillo
INC := -I$(INCDIR) -I/usr/local/include
INCDEP := -I$(INCDIR)

LDFLAGS := `root-config --cflags` `root-config --libs` -pthread

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT)) # get list of all *.cpp files in directory ./src
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT))) # replace PATTERN,REPLACEMENT,TEXT
# the : = operator replaces the file extension ".cpp" with ".o" (just replaces text)
# so we have a list of all *.o files in the format ./src/*.o (instead of ./src/*.cpp
# the patsubst (pattern substitution) replaces the "src" with "obj"
# so SOURCES becomes a list of ./src/*.cpp
# and OBJECTS becomes a list of ./obj/*.o
# 1-to-1 correspondance between objects and source files


#Defauilt Make
#all: resources $(TARGET)
all: $(TARGET)

#Remake
remake: cleaner all

#Copy Resources from Resources Directory to Target Directory
#resources: directories
#	@cp $(RESDIR)/* $(TARGETDIR)/

#Make the Directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)

#Clean only Objecst
clean:
	@$(RM) -rf $(BUILDDIR)

#Full Clean, Objects and Binaries
cleaner: clean
#	@$(RM) -rf $(TARGETDIR)

#Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)

#Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
	@$(CC) $(CFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp

default: build
	@echo "default"

