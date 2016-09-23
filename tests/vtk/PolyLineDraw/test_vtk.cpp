#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
 
int main(int, char *[])
{
    // Create five points. 
    double origin[3] = {0.0, 0.0, 0.0};
    double p0[3] = {1.0, 0.0, 0.0};
    double p1[3] = {0.0, 1.0, 0.0};
    double p2[3] = {0.0, 1.0, 2.0};
    double p3[3] = {1.0, 2.0, 3.0};
 
    // Create a vtkPoints object and store the points in it
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(origin);
    points->InsertNextPoint(p0);
    points->InsertNextPoint(p1);
    points->InsertNextPoint(p2);
    points->InsertNextPoint(p3);
 
    // Create a cell array to store the lines in and add the lines to it
    vtkSmartPointer<vtkCellArray> lines =
        vtkSmartPointer<vtkCellArray>::New();
 
    for(unsigned int i = 0; i < 3; i++)
    {
        vtkSmartPointer<vtkLine> line =
            vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0,i);
        line->GetPointIds()->SetId(1,i+1);
        lines->InsertNextCell(line);
    }
 
    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> linesPolyData =
        vtkSmartPointer<vtkPolyData>::New();
 
    // Add the points to the dataset
    linesPolyData->SetPoints(points);
 
    // Add the lines to the dataset
    linesPolyData->SetLines(lines);
 
    // Setup actor and mapper
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(linesPolyData);
 
    vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
 
    // Setup render window, renderer, and interactor
    vtkSmartPointer<vtkRenderer> renderer =
        vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow =
        vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = 
        vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
    renderer->AddActor(actor);
 
    renderWindow->Render();
    renderWindowInteractor->Start();
 
    return EXIT_SUCCESS;
}
