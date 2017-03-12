/**
 * @file vtk.cpp
 * @brief Contains definitions of classes/functions declared in
 * vtk.hpp .
 */

# include "vtk.hpp"
# include <vtkGenericCell.h>
# include <vtkTubeFilter.h>
# include <vtkLineSource.h>
# include <vtkPolyDataMapper.h>
# include <vtkSphereSource.h>
# include <vtkCylinderSource.h>
# include <vtkAppendFilter.h>
# include <vtkTransform.h>
# include <vtkTransformPolyDataFilter.h>
# include <cmath>


namespace eva {
//------------------------------------------ Contents -------------------------------------------//

vtk_sptr<vtkRenderWindowInteractor>
display(const vtk_sptr<vtkUnstructuredGrid> ugrid, bool deformed);

void
write_vtu(const vtk_sptr<vtkUnstructuredGrid> ugrid,
          const std::string& filename);

void write_3d(const vtk_sptr<vtkUnstructuredGrid> ugrid,
               const std::string& filename);

vtk_sptr<vtkPolyData> to_polydata(vtk_sptr<vtkUnstructuredGrid> ugrid);


//-----------------------------------------------------------------------------------------------//
vtk_sptr<vtkRenderWindowInteractor>
display_grid(const vtk_sptr<vtkUnstructuredGrid> ugrid, bool displace)
{
    // Create a renderer, render window, and interactor
    auto renderer = vtk_sptr<vtkRenderer>::New();
    renderer->SetBackground(0,0,0);

    // Recover cross section information
    auto cell_data = ugrid->GetCellData();
    auto sections  = cell_data->GetScalars("A [m^2]");
    
    // double range[2] = {1., 1.};
    // sections->GetRange(range);
    // auto max_section = range[1];
    constexpr auto alpha = 10.0;
    
    // Recover displacements
    auto point_data = ugrid->GetPointData();
    auto displacements = point_data->GetVectors("Displacements [m]");
    // cell_data->PrintSelf(std::cout, vtkIndent());
    
    // Iterate over all cells
    vtkCellIterator * cell_it = ugrid->NewCellIterator();
    for (cell_it->InitTraversal(); !cell_it->IsDoneWithTraversal(); cell_it->GoToNextCell())
    {
        // Get cell points & build line source out of them
        auto points = cell_it->GetPoints();
        auto line = vtk_sptr<vtkLineSource>::New();
        // if (displace)
        // {
        //     double p1 [] = {0., 0., 0.};
        //     double p2 [] = {0., 0., 0.};
            
        //     points->GetPoint(0, p1);
        //     points->GetPoint(1, p2);

        //     // auto d1 = displacements->GetTuple2(cell_it->GetCellId());
        // }
        // else // Just add points as they are
        {
            line->SetPoints(points);
        }
        
        // Build LINE mapper & actor
        auto line_mapper = vtk_sptr<vtkPolyDataMapper>::New();
        line_mapper->SetInputConnection(line->GetOutputPort());

        auto line_actor = vtk_sptr<vtkActor>::New();
        line_actor->SetMapper(line_mapper);
        line_actor->GetProperty()->SetLineWidth(1.);

        // Build TUBE filter around the line
        auto tube_filter = vtkSmartPointer<vtkTubeFilter>::New();
        tube_filter->SetInputConnection(line->GetOutputPort());
        tube_filter->SetNumberOfSides(50);
        tube_filter->Update();

        // Set TUBE radius
        auto section  = sections->GetTuple1(cell_it->GetCellId());
        auto tube_rad = section * alpha; /// max_section * alpha;  
        tube_filter->SetRadius(tube_rad); //default is .5
        
        // Build TUBE mapper & actor and attach 'em to the line
        auto tube_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        tube_mapper->SetInputConnection(tube_filter->GetOutputPort());
        auto tube_actor = vtkSmartPointer<vtkActor>::New();
        // tube_actor->GetProperty()->SetOpacity(0.5); //Make the tube have some transparency.
        tube_actor->SetMapper(tube_mapper);

        // Build SPHEREs
        std::vector<vtkSmartPointer<vtkActor> > sphere_actors;
 
        for(unsigned int i = 0; i < 2; i++)
        {
            // Build source
            auto sphere = vtkSmartPointer<vtkSphereSource>::New();
            sphere->SetCenter(points->GetPoint(i));
            sphere->SetRadius(2*tube_rad);
            sphere->SetThetaResolution(20);
            sphere->SetPhiResolution(20);
            
            // Build SPHERE mapper & actor
            auto sphere_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            sphere_mapper->SetInputConnection(sphere->GetOutputPort());
            auto sphere_actor = vtkSmartPointer<vtkActor>::New();
            sphere_actor->SetMapper(sphere_mapper);
            
            sphere_actors.push_back(sphere_actor);
        }
        
        // Add actors to renderer
        renderer->AddActor(tube_actor);
        renderer->AddActor(line_actor);
        for (auto & sphere_actor : sphere_actors)
            renderer->AddActor(sphere_actor);
    }

    // Attach window & window interactor
    auto window = vtk_sptr<vtkRenderWindow>::New();
    window->AddRenderer(renderer);

    auto win_interactor = vtk_sptr<vtkRenderWindowInteractor>::New();
    win_interactor->SetRenderWindow(window);    
 
    // Render and interact
    window->Render();
    win_interactor->Start();

    return win_interactor;
}


// vtk_sptr<vtkPolyData> to_polydata(vtk_sptr<vtkUnstructuredGrid> ugrid)
// {
//     auto surface_filter = vtk_sptr<vtkDataSetSurfaceFilter>::New();

//     surface_filter->SetInputData(ugrid);
//     surface_filter->Update(); 
 
//     return surfaceFilter->GetOutput();
// }


void write_3d(const vtk_sptr<vtkUnstructuredGrid> ugrid,
               const std::string& filename)
{
    // Setup stuff to write
    auto polydatas = std::vector< vtk_sptr<vtkPolyData> >();
    
    // Recover cross section information
    auto cell_data = ugrid->GetCellData();
    auto sections  = cell_data->GetScalars("A [m^2]");
    
    double range[2] = {1., 1.};
    sections->GetRange(range);
    auto max_section = range[1];
    constexpr auto alpha = 0.1;
    
    // Recover displacements
    auto point_data = ugrid->GetPointData();
    auto displacements = point_data->GetVectors("Displacements [m]");
    // cell_data->PrintSelf(std::cout, vtkIndent());
    
    // Iterate over all cells
    auto cell_it = ugrid->NewCellIterator();
    for (cell_it->InitTraversal(); !cell_it->IsDoneWithTraversal(); cell_it->GoToNextCell())
    {
        // Build LINE
        auto points = cell_it->GetPoints();
        auto line = vtk_sptr<vtkLineSource>::New();
        line->SetPoints(points);
        // polydatas.push_back(line->GetOutput());

        // Build CYLINDER
        auto cylinder = vtk_sptr<vtkCylinderSource>::New();

        // Set cylinder radius
        auto section  = sections->GetTuple1(cell_it->GetCellId());
        auto cyl_rad = section / max_section * alpha;  
        cylinder->SetRadius(cyl_rad); //default is .5

        // Set cylinder height
        fixed_vector<3> p1(points->GetPoint(0));
        fixed_vector<3> p2(points->GetPoint(1));
        auto height = (p1 - p2).norm();
        cylinder->SetHeight(height);

        cylinder->SetResolution(100);
        cylinder->Update();

        // Rotate & translate cylinder
        auto transform = vtkSmartPointer<vtkTransform>::New();
        transform->PostMultiply();

        fixed_vector<3> dir = (p2 - p1)/(p2 - p1).norm();
        fixed_vector<3> e_y = {0., 1., 0.};
        double angle_deg  = acos(dir.dot(e_y)) * 180 / M_PI;
        transform->RotateZ(angle_deg);

        fixed_vector<3> mid_point = (p1 + p2) / 2.;
        transform->Translate(mid_point.data());
        
        auto transformation = vtk_sptr<vtkTransformPolyDataFilter>::New();
        transformation->SetInputData(cylinder->GetOutput());
        transformation->SetTransform(transform);
        transformation->Update();
        
        polydatas.push_back(transformation->GetOutput());
        
        for(auto i = 0u; i < 2u; ++i)
        {
            // Build sphere
            auto sphere = vtk_sptr<vtkSphereSource>::New();
            
            sphere->SetCenter(points->GetPoint(i));
            sphere->SetRadius(4.);
            sphere->SetThetaResolution(20);
            sphere->SetPhiResolution(20);
            
            polydatas.push_back(sphere->GetOutput());
        }
    }
    
    // Combine datasets
    auto append_filter = vtk_sptr<vtkAppendFilter>::New();
    
    append_filter->AddInputData(ugrid);
    for (auto & polydata : polydatas)
        append_filter->AddInputData(polydata); 

    append_filter->Update();
 
    vtkSmartPointer<vtkUnstructuredGrid> combined = append_filter->GetOutput();
    combined->GetPointData()->AddArray(displacements);

    // Attach data to writer
    auto writer= vtk_sptr<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(combined);
    writer->Write();
}


// vtk_sptr<vtkRenderWindowInteractor>
// display(const vtk_sptr<vtkUnstructuredGrid> ugrid)
// {
//     // Build Mapper
//     vtk_sptr<vtkDataSetMapper> mapper = vtk_sptr<vtkDataSetMapper>::New();
//     mapper->SetInputData(ugrid);

//     // Connect Actor (adjusts visible properties)
//     vtk_sptr<vtkActor> actor = vtk_sptr<vtkActor>::New();
//     actor->SetMapper(mapper);
//     actor->GetProperty()->SetLineWidth(4.);

//     // Connect Renderer and Window
//     vtk_sptr<vtkRenderer> renderer = vtk_sptr<vtkRenderer>::New();
//     vtk_sptr<vtkRenderWindow> renderWindow = vtk_sptr<vtkRenderWindow>::New();
//     renderWindow->AddRenderer(renderer);
//     vtk_sptr<vtkRenderWindowInteractor> renderWindowInteractor =
//         vtk_sptr<vtkRenderWindowInteractor>::New();
//     renderWindowInteractor->SetRenderWindow(renderWindow);

//     // Connect renderer to the actors
//     // for(auto& actor : actors) renderer->AddActor(actor);
//     renderer->AddActor(actor);
//     renderer->SetBackground(0., 0., 0.);
 
//     renderWindow->Render();
//     renderWindowInteractor->Start();

//     return renderWindowInteractor;
// }


void write_vtu(const vtk_sptr<vtkUnstructuredGrid> ugrid,
               const std::string& filename)
{
    auto writer = vtk_sptr<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(ugrid);
    writer->Write();
}


} // end namespace
