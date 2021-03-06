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
# include <vtkXMLPolyDataWriter.h>
# include <vtkUnstructuredGrid.h>
# include <vtkCellIterator.h>
# include <vtkPointData.h>
# include <vtkProperty.h>
# include <vtkCellData.h>
# include <vtkDoubleArray.h>
# include <vtkIndent.h>



namespace eva {

//####################################### DECLARATIONS #############################################

template <typename T> using vtk_sptr = vtkSmartPointer<T>;

template <typename S>
vtk_sptr<vtkUnstructuredGrid> to_unstructured_grid(const S& s);

vtk_sptr<vtkPolyData> to_polydata(vtk_sptr<vtkUnstructuredGrid> ugrid);


template <typename Structure>
vtk_sptr<vtkRenderWindowInteractor>
display(const Structure& structure);

template <typename Structure>
vtk_sptr<vtkRenderWindowInteractor>
display(
    const Structure& structure,
    const std::vector<typename result_of<Structure>::type>& results);

vtk_sptr<vtkRenderWindowInteractor>
display_grid(const vtk_sptr<vtkUnstructuredGrid> ugrid, bool displace = false);


template<typename Structure>
void
write_vtu(const Structure& structure,
          const std::vector<typename result_of<Structure>::type>& results,
          const std::string& filename);

void
write_vtu(const vtk_sptr<vtkUnstructuredGrid> ugrid,
          const std::string& filename);


void write_3d(const vtk_sptr<vtkUnstructuredGrid> ugrid,
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


/// Adds the joint (node) properties and results of the structure to the vtk
/// unstructured grid object
// template <typename S, typename R> 
// void vtk_add_joint_properties(
//     const S& structure,
//     const std::vector<R>& results,
//     vtk_sptr<vtkUnstructuredGrid> ugrid);

/// Adds results (i.e. computed nodal values) to the given unstructured grid
template <typename R>
void vtk_add_joint_results(
    const std::vector<R> & results,
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

#ifdef __EVA_THERMO__

/// Specialization for a truss_element
template <typename BaseKind>
struct vtk_joint_properties_adder< thermo_joint<BaseKind> >;

/// Specialization for a truss_element
template <typename BaseKind>
struct vtk_element_properties_adder< thermo_element<BaseKind> >;


/// Specialization for a truss result type
template <typename BaseKind>
struct vtk_joint_properties_adder< result< thermo_kind<BaseKind> > >;

#endif//__EVA_THERMO__

//####################################### DEFINITIONS ##############################################


template <typename S>
vtk_sptr<vtkUnstructuredGrid> to_vtk_unstructured_grid(
    const S& structure,
    const std::vector<typename result_of<S>::type>& results)
{
    // Build grid and add data
    auto ugrid = to_vtk_unstructured_grid(structure);
    
    // Add joint (node) results
    vtk_add_joint_results(results, ugrid);
    
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

    // Add joint (node) properties
    vtk_add_joint_properties(structure, ugrid);

    // Add element (edge) properties
    vtk_add_element_properties(structure, ugrid);
    
    return ugrid;
}

// ---------- DISPLAY ---------- //

template <typename Structure>
vtk_sptr<vtkRenderWindowInteractor>
display(const Structure& structure)
{
    auto ugrid = to_vtk_unstructured_grid(structure);
    return display_grid(ugrid);
}

template <typename Structure>
vtk_sptr<vtkRenderWindowInteractor>
display(
    const Structure& structure,
    const std::vector<typename result_of<Structure>::type>& results)
{
    auto ugrid = to_vtk_unstructured_grid(structure, results);
    return display_grid(ugrid, true);
}


// ---------- Write VTU ---------- //

template<typename S>
void
write_vtu(const S& structure,
          const std::vector<typename result_of<S>::type>& results,
          const std::string& filename)
{
    auto ugrid = to_vtk_unstructured_grid(structure, results);
    // auto ugrid = to_to_vtk_unstructured_grid(structure);
    
    write_vtu(ugrid, filename);
}

template<typename S>
void
write_vtu(const S& structure, const std::string& filename)
{
    auto ugrid = to_vtk_unstructured_grid(structure);
    write_vtu(ugrid, filename);
}


// ---------- Write 3D ---------- //

template<typename S>
void
write_3d(const S& structure,
          const std::vector<typename result_of<S>::type>& results,
          const std::string& filename)
{
    auto ugrid = to_vtk_unstructured_grid(structure, results);
    write_3d(ugrid, filename);
}

template<typename S>
void
write_3d(const S& structure, const std::string& filename)
{
    auto ugrid = to_vtk_unstructured_grid(structure);
    write_3d(ugrid, filename);
}


// ---------- ADDERS ---------- //

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

// template <typename S, typename R> 
// void vtk_add_joint_properties(
//     const S & structure,
//     const std::vector< typename result_of<S>::type >& results,
//     vtk_sptr<vtkUnstructuredGrid> ugrid)
// {
//     vtk_joint_properties_adder<typename joint_of<S>::type>()(structure, ugrid);
//     // vtk_joint_properties_adder<typename result_of<S>::type>()(structure, results, ugrid);
// }

// template <typename S> 
// void vtk_add_element_properties(
//     const S& structure,
//     const std::vector< typename result_of<S>::type >& results,
//     vtk_sptr<vtkUnstructuredGrid> ugrid)
// {
//     vtk_element_properties_adder<typename element_of<S>::type>()(structure, ugrid);
//     // vtk_element_properties_adder<typename result_of<S>::type>()(structure, results, ugrid);
// }

// template <typename S> 
// void add_joint_properties(
//     const S& structure,
//     const std::vector< typename result_of<S>::type >& results,
//     vtk_sptr<vtkUnstructuredGrid> ugrid)
// {
//     // Add data
//     using vertex_t = typename joint_of<S>::type;
//     vtk_joint_properties_adder<vertex_t>()(structure, ugrid);

//     //Add results
//     using result_t = result<typename kind_of<S>::type>;
//     vtk_joint_properties_adder<result_t>()(structure, results, ugrid);
    
// }


template <typename R>
void vtk_add_joint_results(
    const std::vector<R>& results,
    vtk_sptr<vtkUnstructuredGrid> ugrid)
{
    using vertex_t = R;
    vtk_joint_properties_adder<vertex_t>()(results, ugrid);
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
        // ugrid->GetCellData()->SetScalars(d_A);
    }
};


template <int N>
struct vtk_joint_properties_adder< result< truss_kind<N> > > 
{
    using kind_type = truss_kind<N>;
    
    template <typename R> 
    void operator()(//const S& structure,
                    const std::vector<R>& results,
                    vtk_sptr<vtkUnstructuredGrid> ugrid)
    {        
        // Init results containers
        vtk_sptr<vtkDoubleArray>
            d_displ(vtk_sptr<vtkDoubleArray>::New()), // Displacements
            d_react(vtk_sptr<vtkDoubleArray>::New()); // Reactions

        // Setup names and content dim
        d_displ->SetName("Displacements [m]");
        d_react->SetName("Reactions [N]");
        d_displ->SetNumberOfComponents(kind_type::sdim);
        d_react->SetNumberOfComponents(kind_type::sdim);
        
        // for (const auto& v :
        // boost::make_iterator_range(vertices(structure)))
        for (const auto & vr : results)
        {
            // const auto& vr = results[v];
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
    using kind_type = frame_kind<2>;
    
    template <typename R> 
    void operator()(//const S& structure,
                    const std::vector<R>& results,
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
        d_rotation->SetNumberOfComponents(kind_type::rdim);
        d_torque  ->SetNumberOfComponents(kind_type::rdim);
        
        // for (const auto& v :
        // boost::make_iterator_range(vertices(structure)))
        for (const auto & vr : results)
        { 
            // const auto& vr = results[v];
            d_rotation->InsertNextTupleValue(vr.rotation.data());
            d_torque  ->InsertNextTupleValue(vr.react_torque.data());
        }
        // Add results to the grid
        vtk_joint_properties_adder< result< truss_kind<2> > >()(results, ugrid);
        ugrid->GetPointData()->AddArray(d_rotation);
        ugrid->GetPointData()->AddArray(d_torque);
    }
};
#endif//__EVA_FRAME__


#ifdef __EVA_THERMO__

template<typename BaseKind>
struct vtk_joint_properties_adder< thermo_joint<BaseKind> > 
{    
    template <typename S> 
    void operator()(const S & s, vtk_sptr<vtkUnstructuredGrid> ugrid)
    {
        // Add base kind stuff
        using base_joint_t = typename BaseKind::joint_type;
        vtk_joint_properties_adder<base_joint_t>()(s, ugrid);
        
        // Init results containers
        vtk_sptr<vtkDoubleArray>
            d_T_bc   (vtk_sptr<vtkDoubleArray>::New()),
            d_flux_bc(vtk_sptr<vtkDoubleArray>::New());

        d_T_bc   ->SetName("BC Temperature [K]");
        d_flux_bc->SetName("BC Heat Flux [W/(m*K)]");
        
        for (const auto & v : boost::make_iterator_range(vertices(s)))
        {
            const auto & vp = s[v];
            d_T_bc   ->InsertNextValue(vp.T_bc   );
            d_flux_bc->InsertNextValue(vp.flux_bc);
        }        
        
        // Add properties to grid
        ugrid->GetPointData()->AddArray(d_T_bc   );
        ugrid->GetPointData()->AddArray(d_flux_bc);
    }
};

template<typename BaseKind>
struct vtk_element_properties_adder<thermo_element<BaseKind> > 
{    
    template <typename S> 
    void operator()(const S& s, vtk_sptr<vtkUnstructuredGrid> ugrid)
    {
        // Add base kind stuff
        using base_edge_t = typename BaseKind::element_type;
        vtk_element_properties_adder<base_edge_t>()(s, ugrid);
        
        // Init properties containers
        vtk_sptr<vtkDoubleArray> d_k(vtk_sptr<vtkDoubleArray>::New());
        
        d_k->SetName("Thermal Conductivity [W/(m*K)]");
        
        for (const auto & e : make_iterator_range(edges(s)))
        {
            const auto & ep = s[e];
            d_k->InsertNextValue(ep.k);
        }
        // Add properties to grid
        ugrid->GetCellData()->AddArray(d_k);
    }
};

template<typename BaseKind>
struct vtk_joint_properties_adder< result<thermo_kind<BaseKind> > > 
{
    // using kind_type = thermo_kind<BaseKind>;
    
    template <typename R> 
    void operator()(//const S& structure,
                    const std::vector<R>& results,
                    vtk_sptr<vtkUnstructuredGrid> ugrid)
    {
        // // Add base kind stuff
        // using base_result_t = typename BaseKind::result_type;
        // vtk_joint_properties_adder<base_result_t>()(result, ugrid);
        
        // Init results containers
        vtk_sptr<vtkDoubleArray>
            d_T   (vtk_sptr<vtkDoubleArray>::New()),
            d_flux(vtk_sptr<vtkDoubleArray>::New());

        d_T   ->SetName("Temperature [K]");
        d_flux->SetName("Heat Flux [W/(m*K)]");
        
        // for (const auto& v :
        // boost::make_iterator_range(vertices(structure)))
        for (const auto & rv : results)
        {
            // const auto& rv = results[v];
            d_T   ->InsertNextValue( rv.T  );
            d_flux->InsertNextValue(rv.flux);
        }        
        
        // Add properties to grid
        ugrid->GetPointData()->AddArray(d_T);
        ugrid->GetPointData()->AddArray(d_flux);
    }
};
#endif//__EVA_THERMO__

} // end namespace eva


# endif //__EVA_VTK__
