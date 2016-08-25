# Compiler Setup
CXX      = clang++
CXXFLAGS = -std=c++14 -Wall -g #-O3

# Exacutable name
EXEC = eva

# Include Directories and Files
DIRNAMES = include \
			/usr/include/vtk-6.1
# /usr/include/vtk-6.1/alglib/ /usr/include/vtk-6.1/VPIC/ \
# 			/usr/include/vtk-6.1/vtklibproj4/ /usr/include/vtk-6.1/vtkmetaio/ \
# 			/usr/include/vtk-6.1/vtksys/ /usr/include/vtk-6.1/vtkverdict/
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
LIBPATHS = 	#/usr/lib/paraview \
			/usr/lib/x86_64-linux-gnu/hdf5/serial/lib

LLIBPATHS  = $(foreach dir, $(LIBPATHS), -L$(dir))
RLIBPATHS  = $(foreach dir, $(LIBPATHS), -Wl,-rpath=$(dir))


LIBNAMES = boost_graph
# LIBNAMES = boost_graph vtkIOImport vtkCommonCore vtksys vtkFiltersSources vtkCommonComputationalGeometry vtkCommonDataModel vtkCommonMath vtkCommonMisc vtkCommonSystem vtkCommonTransforms vtkFiltersGeneral vtkFiltersCore vtkCommonExecutionModel vtkRenderingCore vtkFiltersExtraction vtkFiltersStatistics  vtkImagingFourier vtkImagingCore vtkalglib vtkFiltersGeometry hdf5  pthread  z  dl  m  hdf5_hl  jsoncpp  expat  vtkFiltersImaging vtkImagingGeneral vtkImagingSources python2.7  vtkRenderingLOD vtkFiltersModeling vtkIOEnSight vtkFiltersHybrid vtkIOGeometry vtkIOCore vtkIOImage vtkDICOMParser vtkmetaio jpeg  png  tiff  netcdf_c++  netcdf vtkIOXML vtkIOXMLParser vtkImagingHybrid vtkInfovisCore vtkInteractionStyle vtkInteractionWidgets vtkRenderingAnnotation vtkImagingColor vtkRenderingFreeType freetype  vtkftgl vtkRenderingVolume vtkRenderingOpenGL vtkViewsCore vtkFiltersParallelFlowPaths vtkFiltersAMR vtkParallelCore vtkIOLegacy vtkFiltersFlowPaths vtkParallelMPI vtkRenderingGL2PS vtkRenderingContext2D gl2ps  vtkInteractionImage vtkPythonInterpreter sqlite3 vtkTestingRendering xml2  vtkFiltersVerdict vtkIONetCDF vtkRenderingVolumeOpenGL vtkFiltersParallelImaging vtkFiltersParallel Xdmf vtkIOAMR vtkFiltersGeneric vtkGUISupportQt vtkIOPLY vtkCommonColor vtkFiltersTexture vtkFiltersParallelStatistics vtkIOVPIC vtkIOParallelExodus vtkIOExodus vtkexoIIc vtkIOFFMPEG vtkIOMovie theoraenc  theoradec  ogg  vtkRenderingVolumeAMR vtkDomainsChemistry vtkFiltersProgrammable vtkChartsCore vtkRenderingLabel vtkIOLSDyna vtkIOExport vtkIOInfovis vtkRenderingFreeTypeOpenGL vtkFiltersParallelMPI vtkRenderingParallelLIC vtkRenderingLIC vtkIOParallelNetCDF vtkIOParallelLSDyna vtkFiltersHyperTree vtkIOXdmf2 vtkRenderingMatplotlib vtkIOParallel vtkIOMPIImage vtkViewsContext2D vtkImagingMorphological vtkRenderingParallel
LIBS = $(foreach lib, $(LIBNAMES), -l$(lib))


all: $(EXEC)

$(EXEC): $(OBJS)
	@echo " Linking..."
	@mkdir -p bin
	@echo " $(CXX) $^ -o bin/$(EXEC) $(LIBS)"; $(CXX) $(CXXFLAGS) $^ -o bin/$(EXEC) $(LIBS)


$(OBJDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(OBJDIR)
	@echo " $(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<"; $(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<


clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(EXEC)"; $(RM) -r $(EXEC)

distclean:
	@echo " (Dist)Cleaning..."; 
	@echo " $(RM) -r $(OBJDIR) $(EXEC)"; $(RM) -r $(OBJDIR) $(EXEC)
# 	@echo " $(RM) -r $(DBGDIR) $(DBGEXEC)"; $(RM) -r $(DBGDIR) $(DBGEXEC)

parser:
	$(CXX) $(CXXFLAGS) tests/test_parser.cpp -o bin/parser $(INCDIRS)


$(OBJDIR)/test_vtk.o: tests/test_vtk.cpp
	@mkdir -p $(OBJDIR)
	@echo " $(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<"; $(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

vtk: $(OBJDIR)/test_vtk.o
	@echo " Linking..."
	@mkdir -p bin
	@echo " $(CXX) $^ -o bin/test_vtk $(LIBS) $(LLIBPATHS) $(RLIBPATHS)"; $(CXX) $^ -o bin/test_vtk \
												$(LIBS) $(LLIBPATHS) $(RLIBPATHS) $(CXXFLAGS)

eigen:
	$(CXX) $(CXXFLAGS) -o bin/test_eigen tests/eigen.cpp

.PHONY: clean distclean test_armadillo
