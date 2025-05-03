// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
//  -   Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//  -   Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
//      in the documentation and/or other materials provided with the distribution.
//  -   All advertising materials mentioning features or use of this software must display the following acknowledgement:
//          This product includes Kratos Multi-Physics technology.
//  -   Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// #define JSON_INCLUDED

// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/io.h"

#include "includes/model_part_io.h"
#include "includes/reorder_consecutive_model_part_io.h"
#include "includes/gid_io.h"
#include "python/add_io_to_python.h"

#ifdef JSON_INCLUDED
#include "includes/json_io.h"
#endif

namespace Kratos
{

namespace Python
{

typedef ModelPart::DataType DataType;

void WriteNodeMesh( GidIO<>& dummy, GidIO<>::MeshType& rThisMesh )
{
    std::cout<<"start printing nodes mesh "<<std::endl;
    dummy.WriteNodeMesh( rThisMesh );
    std::cout<<"end printing nodes mesh "<<std::endl;
}

void WriteSphereMesh( GidIO<>& dummy, GidIO<>::MeshType& rThisMesh )
{
    //KRATOS_WATCH("writing Sphere Mesh"); //should be a KRATOS_LOG not a KRATOS_WATCH
    dummy.WriteSphereMesh( rThisMesh );
}

void WriteCircleMesh( GidIO<>& dummy, GidIO<>::MeshType& rThisMesh )
{
    //KRATOS_WATCH("writing circle Mesh"); //should be a KRATOS_LOG not a KRATOS_WATCH
    dummy.WriteCircleMesh( rThisMesh );
}

void WriteMesh( GidIO<>& dummy, GidIO<>::MeshType& rThisMesh )
{
    dummy.WriteMesh( rThisMesh );
}

void IntPrintOnGaussPoints( GidIO<>& dummy, const Variable<int>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void DoublePrintOnGaussPoints( GidIO<>& dummy, const Variable<DataType>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void Array1DPrintOnGaussPoints( GidIO<>& dummy, const Variable<array_1d<DataType, 3> >& rVariable,
                                ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void VectorPrintOnGaussPoints( GidIO<>& dummy, const Variable<Vector>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void MatrixPrintOnGaussPoints( GidIO<>& dummy, const Variable<Matrix>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void (GidIO<>::*pointer_to_bool_write_nodal_results)( Variable<bool> const& rVariable,
        GidIO<>::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GidIO<>::WriteNodalResults;
void (GidIO<>::*pointer_to_double_write_nodal_results)( Variable<DataType> const& rVariable,
        GidIO<>::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GidIO<>::WriteNodalResults;
void (GidIO<>::*pointer_to_array1d_write_nodal_results)(
    Variable<array_1d<DataType, 3> > const& rVariable, GidIO<>::NodesContainerType& rNodes,
    double SolutionTag, std::size_t SolutionStepNumber) = &GidIO<>::WriteNodalResults;

//         void (GidIO::*pointer_to_vector_write_nodal_results)(Variable<Vector > const& rVariable, GidIO::NodesContainerType& rNodes, double SolutionTag, std::size_t SolutionStepNumber) = &GidIO::WriteNodalResults;
void (GidIO<>::*pointer_to_matrix_write_nodal_results)(Variable<Matrix > const& rVariable, GidIO<>::NodesContainerType& rNodes, double SolutionTag, std::size_t SolutionStepNumber) = &GidIO<>::WriteNodalResults;

void (GidIO<>::*local_axes_write_nodal_results)( Variable<array_1d<DataType, 3> > const& rVariable,
        GidIO<>::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GidIO<>::WriteLocalAxesOnNodes;

/////////////////////////////////////////////////////////////
/// NON-HISTORICAL DATABASE                               ///
////////////////////////////////////////////////////////////
void (GidIO<>::*pointer_to_bool_write_nodal_results_NH)( Variable<bool> const& rVariable, GidIO<>::NodesContainerType& rNodes, double SolutionTag)
    = &GidIO<>::WriteNodalResultsNonHistorical;

void (GidIO<>::*pointer_to_double_write_nodal_results_NH)( Variable<DataType> const& rVariable,  GidIO<>::NodesContainerType& rNodes, double SolutionTag)
    = &GidIO<>::WriteNodalResultsNonHistorical;
void (GidIO<>::*pointer_to_array1d_write_nodal_results_NH)(Variable<array_1d<DataType, 3> > const& rVariable, GidIO<>::NodesContainerType& rNodes, double SolutionTag)
    = &GidIO<>::WriteNodalResultsNonHistorical;

void (GidIO<>::*pointer_to_matrix_write_nodal_results_NH)(Variable<Matrix > const& rVariable, GidIO<>::NodesContainerType& rNodes, double SolutionTag)
    = &GidIO<>::WriteNodalResultsNonHistorical;

void (GidIO<>::*local_axes_write_nodal_results_NH)( Variable<array_1d<DataType, 3> > const& rVariable, GidIO<>::NodesContainerType& rNodes, double SolutionTag)
    = &GidIO<>::WriteLocalAxesOnNodesNonHistorical;


//         void (GidIO::*pointer_to_double_cond_print_on_gauss_points)(const Variable<double>& rVariable,
//               ModelPart& r_model_part, double SolutionTag) = &GidIO::CondPrintOnGaussPoints;
//         void (GidIO<>::*pointer_to_double_print_on_gauss_points)(const Variable<double >& rVariable,
//               ModelPart& r_model_part, double SolutionTag)
//                 = &GidIO<>::PrintOnGaussPoints;
//         void (GidIO<>::*pointer_to_matrix_print_on_gauss_points)(const Variable<Matrix >& rVariable,
//               ModelPart& r_model_part, double SolutionTag)
//                 = &GidIO<>::PrintOnGaussPoints;
//         void (GidIO<>::*pointer_to_vector_print_on_gauss_points)(const Variable<Vector >& rVariable,
//               ModelPart& r_model_part, double SolutionTag)
//                 = &GidIO<>::PrintOnGaussPoints;

template<class TIOType>
void ReadInitialValues1(TIOType& IO, typename TIOType::NodesContainerType& rThisNodes,
        typename TIOType::ElementsContainerType& rThisElements, typename TIOType::ConditionsContainerType& rThisConditions)
{
    IO.ReadInitialValues(rThisNodes, rThisElements, rThisConditions);
}

template<class TIOType>
void ReadInitialValues2(TIOType& IO, typename TIOType::ModelPartType& rThisModelPart)
{
    IO.ReadInitialValues(rThisModelPart);
}

template<class TModelPartType>
void AddModelPartIOToPython(const std::string& Prefix)
{
    using namespace boost::python;

    typedef IO<TModelPartType> IOType;

    void (IOType::*pointer_to_io_read_single_properties)(Properties& rThisProperties ) = &IOType::ReadProperties;

    void (IOType::*pointer_to_io_read_properties)(typename IOType::PropertiesContainerType& rThisProperties ) = &IOType::ReadProperties;

    auto io_python_interface = class_<IOType, typename IOType::Pointer, bases<BaseIO>, boost::noncopyable>((Prefix+"IO").c_str())
    .def("ReadNode", &IOType::ReadNode)
    .def("ReadNodes", &IOType::ReadNodes)
    .def("WriteNodes", &IOType::WriteNodes)
    .def("ReadProperties", pointer_to_io_read_single_properties)
    .def("ReadProperties", pointer_to_io_read_properties)
    .def("ReadElements", &IOType::ReadElements)
    .def("WriteElements", &IOType::WriteElements)
    .def("ReadConditions", &IOType::ReadConditions)
    .def("ReadInitialValues", &ReadInitialValues1<IOType>)
    .def("ReadInitialValues", &ReadInitialValues2<IOType>)
    .def("ReadMesh", &IOType::ReadMesh)
    .def("ReadModelPart", &IOType::ReadModelPart)
    ;
    io_python_interface.attr("READ") = IOType::READ;
    io_python_interface.attr("WRITE") = IOType::WRITE;
    io_python_interface.attr("APPEND") = IOType::APPEND;
    io_python_interface.attr("IGNORE_VARIABLES_ERROR" ) = IOType::IGNORE_VARIABLES_ERROR;

    typedef ModelPartIO<TModelPartType> ModelPartIOType;
    class_<ModelPartIOType, typename ModelPartIOType::Pointer, bases<IOType>,  boost::noncopyable>(
        (Prefix+"ModelPartIO").c_str(), init<std::string const&>())
        .def(init<std::string const&, const Flags>())
    ;

    typedef ReorderConsecutiveModelPartIO<TModelPartType> ReorderConsecutiveModelPartIOType;
    class_<ReorderConsecutiveModelPartIOType, typename ReorderConsecutiveModelPartIOType::Pointer, bases<ModelPartIOType>, boost::noncopyable>(
        (Prefix+"ReorderConsecutiveModelPartIO").c_str(), init<std::string const&>())
        .def(init<std::string const&, const Flags>())
    ;
}

void AddIOToPython()
{
    using namespace boost::python;

    class_<BaseIO, BaseIO::Pointer, boost::noncopyable>("BaseIO", init<>())
    .def(self_ns::str(self));

    AddModelPartIOToPython<ModelPart>("");
    AddModelPartIOToPython<ComplexModelPart>("Complex");
    AddModelPartIOToPython<GComplexModelPart>("GComplex");

#ifdef JSON_INCLUDED
    class_<KratosJsonIO, KratosJsonIO::Pointer, bases<BaseIO>, boost::noncopyable>(
         "JsonIO",init<std::string const&>())
        .def(init<std::string const&, const Flags>())
    ;
#endif

    class_<GidIO<>, GidIO<>::Pointer, bases<BaseIO>, boost::noncopyable>(
        "GidIO",init<std::string const&, GiD_PostMode,
        MultiFileFlag,
        WriteDeformedMeshFlag,
        WriteConditionsFlag>())
    //.def(init<std::string const&>())
    .def("WriteMesh",WriteMesh)
    .def("WriteNodeMesh",WriteNodeMesh)
    .def("WriteSphereMesh",WriteSphereMesh)
    .def("WriteCircleMesh",WriteCircleMesh)

    .def("InitializeMesh",&GidIO<>::InitializeMesh)
    .def("FinalizeMesh",&GidIO<>::FinalizeMesh)

    .def("InitializeResults",&GidIO<>::InitializeResults)
    .def("FinalizeResults",&GidIO<>::FinalizeResults)

    .def("WriteNodalResults",pointer_to_bool_write_nodal_results)
    .def("WriteNodalResults",pointer_to_double_write_nodal_results)
    .def("WriteNodalResults",pointer_to_array1d_write_nodal_results)

//                    .def("WriteNodalResults",pointer_to_vector_write_nodal_results)
    .def("WriteNodalResults",pointer_to_matrix_write_nodal_results)

    .def("WriteLocalAxesOnNodes",local_axes_write_nodal_results)
    // NonHistorical
    .def("WriteNodalResultsNonHistorical",pointer_to_bool_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical",pointer_to_double_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical",pointer_to_array1d_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical",pointer_to_matrix_write_nodal_results_NH)
    .def("WriteLocalAxesOnNodesNonHistorical",local_axes_write_nodal_results_NH)

//                     .def("PrintOnGaussPoints", pointer_to_double_print_on_gauss_points)
    .def("PrintOnGaussPoints", IntPrintOnGaussPoints)
    .def("PrintOnGaussPoints", DoublePrintOnGaussPoints)
    .def("PrintOnGaussPoints", Array1DPrintOnGaussPoints)
    .def("PrintOnGaussPoints", VectorPrintOnGaussPoints)
    .def("PrintOnGaussPoints", MatrixPrintOnGaussPoints)


//                     .def("PrintOnGaussPoints", pointer_to_vector_print_on_gauss_points)
//                     .def("PrintOnGaussPoints", pointer_to_matrix_print_on_gauss_points)
//                     .def("CondPrintOnGaussPoints", pointer_to_double_cond_print_on_gauss_points)

    .def("Flush",&GidIO<>::Flush)
    .def("ChangeOutputName",&GidIO<>::ChangeOutputName)
    .def("CloseResultFile",&GidIO<>::CloseResultFile)
    //.def("",&DatafileIO::)
    //.def(self_ns::str(self))
    ;

    enum_<GiD_PostMode>("GiDPostMode")
    .value("GiD_PostAscii", GiD_PostAscii)
    .value("GiD_PostAsciiZipped", GiD_PostAsciiZipped)
    .value("GiD_PostBinary", GiD_PostBinary)
    .value("GiD_PostHDF5", GiD_PostHDF5)
    ;

    enum_<WriteDeformedMeshFlag>("WriteDeformedMeshFlag")
    .value("WriteDeformed", WriteDeformed)
    .value("WriteUndeformed", WriteUndeformed)
    ;

    enum_<WriteConditionsFlag>("WriteConditionsFlag")
    .value("WriteConditions",WriteConditions)
    .value("WriteElementsOnly",WriteElementsOnly)
    .value("WriteConditionsOnly",WriteConditionsOnly)
    ;

    enum_<MultiFileFlag>("MultiFileFlag")
    .value("SingleFile",SingleFile)
    .value("MultipleFiles",MultipleFiles)
    ;
}
}  // namespace Python.

} // Namespace Kratos

