#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkPoints2D.h>
#include <vtkVertexGlyphFilter.h>

# include "truss.hpp"
# include "frame.hpp"
# include "io.hpp"


template <typename T> using vtk_ptr = vtkSmartPointer<T>;

template <typename S>
vtk_ptr<vtkUnstructuredGrid> to_unstructured_grid(
    const S& s
    )
{
    // Determine vtk point type accordingly to the nr of spatial dimensions
    constexpr static int sdim = S::graph_bundled::sdim;
    // using points_t = typename std::conditional<
    //     sdim == 3, vtkPoints, vtkPoints2D
    //     >::type;
    using points_t = vtkPoints;
    
    // Build VTK points
    vtk_ptr<points_t> points = vtk_ptr<points_t>::New();
    for (const auto& v : make_iterator_range(vertices(s))) {
        // Build point from joint coordinates
        const auto coords = s[v].coords.data();
        // points->InsertPoint(v, coords);
        if (sdim == 3)
            points->InsertPoint(v, coords);
        else /*sdim = 2 */
            points->InsertPoint(v, coords[0], coords[1], 0.0);    
    }

    // Build VTK cells (lines representing elements)
    vtk_ptr<vtkCellArray> cells = vtk_ptr<vtkCellArray>::New();
    for (const auto& e : make_iterator_range(edges(s))) {
        // Build line (#local dof = {0,1}, #global dof = {vertices})
        vtk_ptr<vtkLine> line = vtk_ptr<vtkLine>::New();
        line->GetPointIds()->SetId(0, source(e, s));
        line->GetPointIds()->SetId(1, target(e, s));
        // Add line to cells
        cells->InsertNextCell(line);
    }

    // Build the unstructured grid
    vtk_ptr<vtkUnstructuredGrid> ugrid = vtk_ptr<vtkUnstructuredGrid>::New();
    ugrid->SetPoints(points);
    ugrid->SetCells(VTK_LINE, cells);

    // Add properties
    // ...

    return ugrid;
}


template <typename Structure>
vtk_ptr<vtkRenderWindowInteractor>
display(const Structure& structure)
{
    auto ugrid = to_unstructured_grid(structure);
    return display(ugrid);
}



vtk_ptr<vtkRenderWindowInteractor>
display(const vtk_ptr<vtkUnstructuredGrid> ugrid)
{
    // Build Mapper
    vtk_ptr<vtkDataSetMapper> mapper = vtk_ptr<vtkDataSetMapper>::New();
    mapper->SetInputData(ugrid);

    // Connect Actor (adjusts visible properties)
    vtk_ptr<vtkActor> actor = vtk_ptr<vtkActor>::New();
    actor->SetMapper(mapper);

    // Connect Renderer and Window
    vtk_ptr<vtkRenderer> renderer = vtk_ptr<vtkRenderer>::New();
    vtk_ptr<vtkRenderWindow> renderWindow = vtk_ptr<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtk_ptr<vtkRenderWindowInteractor> renderWindowInteractor =
        vtk_ptr<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
 
    renderer->AddActor(actor);
    renderer->SetBackground(.3, .6, .3); // Background color green
 
    renderWindow->Render();
    renderWindowInteractor->Start();

    return renderWindowInteractor;
}

    
int main(int argc, char *argv[])
{
    // Structure typedef
    using structure_type = eva::truss3d;
    
    // Read input file and build structure
    auto filename  = std::string("/home/claudio/eva/input/truss3d.dot");
    auto structure = eva::read_from_graphviz<structure_type>(filename);

    auto ugrid = to_unstructured_grid(structure);
 
    // Write file
    vtk_ptr<vtkXMLUnstructuredGridWriter> writer =
        vtk_ptr<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName("test.vtu");
    writer->SetInputData(ugrid);
    writer->Write();
 
    // Read and display file for verification that it was written correclty
    // vtk_ptr<vtkXMLUnstructuredGridReader> reader =
    //     vtk_ptr<vtkXMLUnstructuredGridReader>::New();
    // reader->SetFileName("test.vtu");
    // reader->Update();
    display(structure);

 
    return EXIT_SUCCESS;
}