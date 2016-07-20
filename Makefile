# Compiler Setup
CXX      = clang++
CXXFLAGS = -std=c++14 -Wall -g #-O3

# Exacutable name
EXEC = eva

# Include Directories and Files
DIRNAMES = include
INCDIRS  = $(foreach dir, $(DIRNAMES), -I$(dir))
INCEXT   = hpp
INCS     = $(shell find $(DIRNAMES) -type f -name *.$(INCEXT))


# Sources and Object Files
SRCDIR = src
OBJDIR = obj
SRCEXT = cpp
SRCS   = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJS   = $(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(SRCS:.$(SRCEXT)=.o))

# External Libraries
LIBNAMES = boost_graph
LIBS     = $(foreach lib, $(LIBNAMES), -l$(lib))


all: $(EXEC)

$(EXEC): $(OBJS)
	@echo " Linking..."
	@mkdir -p bin
	@echo " $(CXX) $^ -o bin/$(EXEC) $(LIBS)"; $(CXX) $^ -o bin/$(EXEC) $(LIBS)
	

$(OBJDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(OBJDIR)
	@echo " $(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<"; $(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<


clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(EXEC)"; $(RM) -r $(EXEC)
		
distclean:
	@echo " (Dist)Cleaning..."; 
	@echo " $(RM) -r $(OBJDIR) $(EXEC)"; $(RM) -r $(OBJDIR) $(EXEC)
#~ 	@echo " $(RM) -r $(DBGDIR) $(DBGEXEC)"; $(RM) -r $(DBGDIR) $(DBGEXEC)

parser:
	$(CXX) $(CXXFLAGS) tests/test_parser.cpp -o bin/parser $(INCDIRS)
#~ test_graphviz:
#~ 	$(CXX) $(CXXFLAGS) tests/test_graphviz.cpp -o tests/graph -lboost_graph
	
.PHONY: clean distclean test_armadillo
