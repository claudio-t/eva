# ifndef __EVA_VTK__
# define __EVA_VTK__


#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkPoints2D.h>
// #include <vtkXMLUnstructuredGridReader.h>
// #include <vtkVertexGlyphFilter.h>

namespace eva {

//####################################### DECLARATIONS #############################################

template <typename T> using vtk_sptr = vtkSmartPointer<T>;

template <typename S>
vtk_sptr<vtkUnstructuredGrid> to_unstructured_grid(const S& s);

template <typename Structure>
vtk_sptr<vtkRenderWindowInteractor>
display(const Structure& structure);

vtk_sptr<vtkRenderWindowInteractor>
display(const vtk_sptr<vtkUnstructuredGrid> ugrid);


//####################################### DEFINITIONS ##############################################

template <typename S>
vtk_sptr<vtkUnstructuredGrid> to_vtk_unstructured_grid(const S& s)
{
    // Determine vtk point type accordingly to the nr of spatial dimensions
    constexpr static int sdim = S::graph_bundled::sdim;
    using points_t = vtkPoints;
    
    // Build VTK points
    vtk_sptr<points_t> points = vtk_sptr<points_t>::New();
    for (const auto& v : make_iterator_range(vertices(s))) {
        // Build point from joint coordinates
        const auto coords = s[v].coords.data();

        if (sdim == 3)
            points->InsertPoint(v, coords);
        else /*sdim = 2 */
            points->InsertPoint(v, coords[0], coords[1], 0.0);    
    }

    // Build VTK cells (lines representing elements)
    vtk_sptr<vtkCellArray> cells = vtk_sptr<vtkCellArray>::New();
    for (const auto& e : make_iterator_range(edges(s))) {
        // Build line (#local dof = {0,1}, #global dof = {vertices})
        vtk_sptr<vtkLine> line = vtk_sptr<vtkLine>::New();
        line->GetPointIds()->SetId(0, source(e, s));
        line->GetPointIds()->SetId(1, target(e, s));
        // Add line to cells
        cells->InsertNextCell(line);
    }

    // Build the unstructured grid
    vtk_sptr<vtkUnstructuredGrid> ugrid = vtk_sptr<vtkUnstructuredGrid>::New();
    ugrid->SetPoints(points);
    ugrid->SetCells(VTK_LINE, cells);

    // Add properties
    // ...

    return ugrid;
}


template <typename Structure>
vtk_sptr<vtkRenderWindowInteractor>
display(const Structure& structure)
{
    auto ugrid = to_vtk_unstructured_grid(structure);
    return display(ugrid);
}



vtk_sptr<vtkRenderWindowInteractor>
display(const vtk_sptr<vtkUnstructuredGrid> ugrid)
{
    // Build Mapper
    vtk_sptr<vtkDataSetMapper> mapper = vtk_sptr<vtkDataSetMapper>::New();
    mapper->SetInputData(ugrid);

    // Connect Actor (adjusts visible properties)
    vtk_sptr<vtkActor> actor = vtk_sptr<vtkActor>::New();
    actor->SetMapper(mapper);

    // Connect Renderer and Window
    vtk_sptr<vtkRenderer> renderer = vtk_sptr<vtkRenderer>::New();
    vtk_sptr<vtkRenderWindow> renderWindow = vtk_sptr<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtk_sptr<vtkRenderWindowInteractor> renderWindowInteractor =
        vtk_sptr<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
 
    renderer->AddActor(actor);
    renderer->SetBackground(.3, .6, .3); // Background color green
 
    renderWindow->Render();
    renderWindowInteractor->Start();

    return renderWindowInteractor;
}



} // end namespace eva


# endif //__EVA_VTK__
