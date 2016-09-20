/**
 * @file vtk.cpp
 * @brief Contains definitions of classes/functions declared in
 * vtk.hpp .
 */

# include "vtk.hpp"

namespace eva {
//------------------------------------------ Contents -------------------------------------------//

vtk_sptr<vtkRenderWindowInteractor>
display(const vtk_sptr<vtkUnstructuredGrid> ugrid);

void
write_vtu(const vtk_sptr<vtkUnstructuredGrid> ugrid,
          const std::string& filename);


//-----------------------------------------------------------------------------------------------//
vtk_sptr<vtkRenderWindowInteractor>
display(const vtk_sptr<vtkUnstructuredGrid> ugrid)
{
    // Build Mapper
    vtk_sptr<vtkDataSetMapper> mapper = vtk_sptr<vtkDataSetMapper>::New();
    mapper->SetInputData(ugrid);

    // Connect Actor (adjusts visible properties)
    vtk_sptr<vtkActor> actor = vtk_sptr<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetLineWidth(4.);

    // Connect Renderer and Window
    vtk_sptr<vtkRenderer> renderer = vtk_sptr<vtkRenderer>::New();
    vtk_sptr<vtkRenderWindow> renderWindow = vtk_sptr<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtk_sptr<vtkRenderWindowInteractor> renderWindowInteractor =
        vtk_sptr<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Connect renderer to the actors
    // for(auto& actor : actors) renderer->AddActor(actor);
    renderer->AddActor(actor);
    renderer->SetBackground(0., 0., 0.);
 
    renderWindow->Render();
    renderWindowInteractor->Start();

    return renderWindowInteractor;
}


void write_vtu(const vtk_sptr<vtkUnstructuredGrid> ugrid,
               const std::string& filename)
{
    vtk_sptr<vtkXMLUnstructuredGridWriter> writer =
        vtk_sptr<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(ugrid);
    writer->Write();
}


} // end namespace
