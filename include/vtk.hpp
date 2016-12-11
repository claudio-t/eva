# ifndef __EVA_VTK__
# define __EVA_VTK__
/**
 * @file vtk.hpp
 * @brief Contains classes/functions used for displaying and saving
 * each supported structure exploiting vtk.
 */

// eva
// # include "truss.hpp"
// # include "frame.hpp"
# include "core.hpp"
// vtk
# include <vtkVersion.h>
# include <vtkSmartPointer.h>
# include <vtkLine.h>
# include <vtkCellArray.h>
# include <vtkDataSetMapper.h>
# include <vtkActor.h>
# include <vtkRenderer.h>
# include <vtkRenderWindow.h>
# include <vtkRenderWindowInteractor.h>
# include <vtkXMLUnstructuredGridWriter.h>
# include <vtkUnstructuredGrid.h>
# include <vtkPointData.h>
# include <vtkProperty.h>
# include <vtkCellData.h>
# include <vtkDoubleArray.h>


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


template<typename S>
void
write_vtu(const S& structure,
          const std::vector<typename result_of<S>::type>& results,
          const std::string& filename);

void
write_vtu(const vtk_sptr<vtkUnstructuredGrid> ugrid,
          const std::string& filename);


/// Adds the joint (node) properties of the structure to the vtk
/// unstructured grid object 
template <typename S>
void vtk_add_joint_properties(const S& structure,
                              vtk_sptr<vtkUnstructuredGrid> ugrid);

/// Adds the element (edge) properties of the structure to the vtk
/// unstructured grid object 
template <typename S>
void vtk_add_element_properties(const S& structure,
                                vtk_sptr<vtkUnstructuredGrid> ugrid);


/// Adds the results of the solved structure to the vtk
/// unstructured grid object 
template <typename S> 
void add_results(
    const S& structure,
    const std::vector< typename result_of<S>::type >& results,
    vtk_sptr<vtkUnstructuredGrid> ugrid);



/// Adds the joint (node) properties and results of the structure to the vtk
/// unstructured grid object
template <typename S> 
void vtk_add_joint_properties(
    const S& structure,
    const std::vector< typename result_of<S>::type >& results,
    vtk_sptr<vtkUnstructuredGrid> ugrid);

/// Adds the joint (node) properties and results of the structure to the vtk
/// unstructured grid object
template <typename S> 
void vtk_add_joint_properties(
    const S& structure,
    const std::vector< typename result_of<S>::type >& results,
    vtk_sptr<vtkUnstructuredGrid> ugrid);


/// Functor for adding joint properties to the vtk unstructured grid object
template <typename JointProps>
struct vtk_joint_properties_adder;

/// Functor for adding element properties to the vtk unstructured grid object
template <typename ElementProps>
struct vtk_element_properties_adder;

#ifdef __EVA_TRUSS__

/// Specialization for a truss_element
template <int N>
struct vtk_joint_properties_adder< truss_joint<N> >;

/// Specialization for a truss_element
template <int N>
struct vtk_element_properties_adder< truss_element<N> >;


/// Specialization for a truss result type
template <int N>
struct vtk_joint_properties_adder< result< truss_kind<N> > >;

#endif//__EVA_TRUSS__


#ifdef __EVA_FRAME__

/// Specialization for a 2D frame result type
template <>
struct vtk_joint_properties_adder< result< frame_kind<2> > >;

/// Specialization for a 2D frame element
template <>
struct vtk_joint_properties_adder< frame_joint<2> >;

/// Specialization for a 2D frame element
template <>
struct vtk_element_properties_adder< frame_element<2> >;

#endif//__EVA_FRAME__

//####################################### DEFINITIONS ##############################################


template <typename S>
vtk_sptr<vtkUnstructuredGrid> to_vtk_unstructured_grid(
    const S& structure,
    const std::vector<typename result_of<S>::type>& results)
{
    // Build grid and add data
    auto ugrid = to_vtk_unstructured_grid(structure);
    
    // Add joint (node) results
    vtk_add_joint_properties(structure, results, ugrid);

    // Add element (edge) properties
    // vtk_add_element_properties(s, ugrid);
    
    return ugrid;
}

template <typename S>
vtk_sptr<vtkUnstructuredGrid> to_vtk_unstructured_grid(const S& structure)
{
    // Determine vtk point type accordingly to the nr of spatial dimensions
    constexpr static int sdim = kind_of<S>::type::sdim;
    using points_t = vtkPoints;
    
    // Build VTK points
    vtk_sptr<points_t> points = vtk_sptr<points_t>::New();
    for (const auto& v : boost::make_iterator_range(vertices(structure)))
    {
        // Build point from joint coordinates
        const auto coords = structure[v].coords.data();
        if (sdim == 3)
            points->InsertPoint(v, coords);
        else /*sdim = 2 */
            points->InsertPoint(v, coords[0], coords[1], 0.0);    
    }

    // Build VTK cells (lines representing elements)
    vtk_sptr<vtkCellArray> cells = vtk_sptr<vtkCellArray>::New();
    for (const auto& e : make_iterator_range(edges(structure)))
    {
        // Build line (#local dof = {0,1}, #global dof = {vertices})
        vtk_sptr<vtkLine> line = vtk_sptr<vtkLine>::New();
        line->GetPointIds()->SetId(0, source(e, structure));
        line->GetPointIds()->SetId(1, target(e, structure));
        
        // Add line to cells
        cells->InsertNextCell(line);
    }

    // Build the unstructured grid
    vtk_sptr<vtkUnstructuredGrid> ugrid = vtk_sptr<vtkUnstructuredGrid>::New();
    ugrid->SetPoints(points);
    ugrid->SetCells(VTK_LINE, cells);

    // Add element (edge) properties
    vtk_add_joint_properties(structure, ugrid);

    // Add element (edge) properties
    vtk_add_element_properties(structure, ugrid);
    
    return ugrid;
}


template <typename Structure>
vtk_sptr<vtkRenderWindowInteractor>
display(const Structure& structure)
{
    auto ugrid = to_vtk_unstructured_grid(structure);
    return display(ugrid);
}


template<typename S>
void
write_vtu(const S& structure,
          const std::vector<typename result_of<S>::type>& results,
          const std::string& filename)
{
    auto ugrid = to_vtk_unstructured_grid(structure, results);
    write_vtu(ugrid, filename);
}

template<typename S>
void
write_vtu(const S& structure, const std::string& filename)
{
    auto ugrid = to_vtk_unstructured_grid(structure);
    write_vtu(ugrid, filename);
}


template <typename S> 
void vtk_add_joint_properties(const S& structure, vtk_sptr<vtkUnstructuredGrid> ugrid)
{
    vtk_joint_properties_adder<typename joint_of<S>::type>()(structure, ugrid);
}

template <typename S> 
void vtk_add_element_properties(const S& structure, vtk_sptr<vtkUnstructuredGrid> ugrid)
{
    vtk_element_properties_adder<typename element_of<S>::type>()(structure, ugrid);
}

template <typename S> 
void vtk_add_joint_properties(
    const S& structure,
    const std::vector< typename result_of<S>::type >& results,
    vtk_sptr<vtkUnstructuredGrid> ugrid)
{
    vtk_joint_properties_adder<typename joint_of<S>::type>()(structure, ugrid);
    vtk_joint_properties_adder<typename result_of<S>::type>()(structure, results, ugrid);
}

template <typename S> 
void vtk_add_element_properties(
    const S& structure,
    const std::vector< typename result_of<S>::type >& results,
    vtk_sptr<vtkUnstructuredGrid> ugrid)
{
    vtk_element_properties_adder<typename element_of<S>::type>()(structure, ugrid);
    // vtk_element_properties_adder<typename result_of<S>::type>()(structure, results, ugrid);
}

template <typename S> 
void add_joint(
    const S& structure,
    const std::vector< typename result_of<S>::type >& results,
    vtk_sptr<vtkUnstructuredGrid> ugrid)
{
    // Add data
    using vertex_t = typename joint_of<S>::type;
    vtk_joint_properties_adder<vertex_t>()(structure, ugrid);

    //Add results
    using result_t = result<typename kind_of<S>::type>;
    vtk_joint_properties_adder<result_t>()(structure, results, ugrid);
    
}

#ifdef __EVA_TRUSS__

template <int N>
struct vtk_joint_properties_adder< truss_joint<N> > 
{    
    template <typename S> 
    void operator()(const S& s, vtk_sptr<vtkUnstructuredGrid> ugrid)
    {        
        // Init properties containers
        vtk_sptr<vtkDoubleArray> d_load  (vtk_sptr<vtkDoubleArray>::New()); // Applied
                                                                            // load
        d_load->SetName("Applied Load [N]");
        d_load->SetNumberOfComponents(kind_of<S>::type::sdim);
        
        for (const auto& v : boost::make_iterator_range(vertices(s)))
        {
            const auto& vp = s[v];
            d_load->InsertNextTupleValue(vp.load.data());
        }
        // Add properties to the grid
        ugrid->GetPointData()->AddArray(d_load);
    }
};

template <int N>
struct vtk_element_properties_adder< truss_element<N> > 
{    
    template <typename S> 
    void operator()(const S& s, vtk_sptr<vtkUnstructuredGrid> ugrid)
    {
        // Init properties containers
        vtk_sptr<vtkDoubleArray>
            d_E(vtk_sptr<vtkDoubleArray>::New()), // Young's modulus
            d_A(vtk_sptr<vtkDoubleArray>::New()); // Cross sectional
                                                  // area
        d_E->SetName("E [Pa]");
        d_A->SetName("A [m^2]");
        
        for (const auto& e : make_iterator_range(edges(s)))
        {
            const auto& ep = s[e];
            d_E->InsertNextValue(ep.E);
            d_A->InsertNextValue(ep.A);
        }
        // Add properties to grid
        ugrid->GetCellData()->AddArray(d_E);
        ugrid->GetCellData()->AddArray(d_A);
    }
};


template <int N>
struct vtk_joint_properties_adder< result< truss_kind<N> > > 
{    
    template <typename S> 
    void operator()(const S& structure,
                    const std::vector< typename result_of<S>::type >& results,
                    vtk_sptr<vtkUnstructuredGrid> ugrid)
    {        
        // Init results containers
        vtk_sptr<vtkDoubleArray>
            d_displ(vtk_sptr<vtkDoubleArray>::New()), // Displacements
            d_react(vtk_sptr<vtkDoubleArray>::New()); // Reactions

        // Setup names and content dim
        d_displ->SetName("Displacements [m]");
        d_react->SetName("Reactions [N]");
        d_displ->SetNumberOfComponents(kind_of<S>::type::sdim);
        d_react->SetNumberOfComponents(kind_of<S>::type::sdim);
        
        for (const auto& v : boost::make_iterator_range(vertices(structure)))
        {
            const auto& vr = results[v];
            d_displ->InsertNextTupleValue(vr.displacement.data());
            d_react->InsertNextTupleValue(vr.reaction.data());
        }
        // Add results to the grid
        ugrid->GetPointData()->AddArray(d_displ);
        ugrid->GetPointData()->AddArray(d_react);
    }
};

#endif//__EVA_TRUSS__


#ifdef __EVA_FRAME__

template <>
struct vtk_joint_properties_adder< frame_joint<2> > 
{    
    template <typename S> 
    void operator()(const S& s, vtk_sptr<vtkUnstructuredGrid> ugrid)
    {   
        // Init properties containers
        vtk_sptr<vtkDoubleArray> d_torque(vtk_sptr<vtkDoubleArray>::New()); // Applied
                                                                            // torque
        d_torque->SetName("Applied Torque [Nm]");
        d_torque->SetNumberOfComponents(kind_of<S>::type::rdim);
        
        for (const auto& v : boost::make_iterator_range(vertices(s)))
        {
            const auto& vp = s[v];
            d_torque->InsertNextTupleValue(vp.torque.data());
        }
        // Add properties to the grid
        vtk_joint_properties_adder<truss_joint<2>>()(s, ugrid);
        ugrid->GetPointData()->AddArray(d_torque);
    }
};

template <>
struct vtk_element_properties_adder< frame_element<2> > 
{    
    template <typename S> 
    void operator()(const S& s, vtk_sptr<vtkUnstructuredGrid> ugrid)
    {   
        // Init properties containers
        vtk_sptr<vtkDoubleArray>
            d_I(vtk_sptr<vtkDoubleArray>::New()); // Second moment of
                                                  // area
        d_I->SetName("I [m^4]");
    
        for (const auto& e : make_iterator_range(edges(s)))
            d_I->InsertNextValue(s[e].I);
        
        // Add properties to grid
        vtk_element_properties_adder< truss_element<2> >()(s, ugrid);
        ugrid->GetCellData()->AddArray(d_I);
    }
};

template <>
struct vtk_joint_properties_adder< result< frame_kind<2> > > 
{    
    template <typename S> 
    void operator()(const S& structure,
                    const std::vector< typename result_of<S>::type >& results,
                    vtk_sptr<vtkUnstructuredGrid> ugrid)
    {
        // Init results containers
        vtk_sptr<vtkDoubleArray>
            d_rotation(vtk_sptr<vtkDoubleArray>::New()), // Rotations
            d_torque  (vtk_sptr<vtkDoubleArray>::New()); // Reaction
                                                         // torques
        // Setup names and content dim
        d_rotation->SetName("Rotations [m]");        
        d_torque  ->SetName("Reaction Torques [Nm]");
        d_rotation->SetNumberOfComponents(kind_of<S>::type::rdim);
        d_torque  ->SetNumberOfComponents(kind_of<S>::type::rdim);
        
        for (const auto& v : boost::make_iterator_range(vertices(structure)))
        {
            const auto& vr = results[v];
            d_rotation->InsertNextTupleValue(vr.rotation.data());
            d_torque  ->InsertNextTupleValue(vr.react_torque.data());
        }
        // Add results to the grid
        vtk_joint_properties_adder< result< truss_kind<2> > >()(structure, results, ugrid);
        ugrid->GetPointData()->AddArray(d_torque);
    }
};

#endif//__EVA_FRAME__

} // end namespace eva


# endif //__EVA_VTK__
