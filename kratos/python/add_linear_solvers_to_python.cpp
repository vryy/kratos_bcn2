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


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "python/add_equation_systems_to_python.h"
#include "linear_solvers/cg_solver.h"
#include "linear_solvers/deflated_cg_solver.h"
#include "linear_solvers/bicgstab_solver.h"
#include "linear_solvers/tfqmr_solver.h"
#include "includes/dof.h"
#include "spaces/ublas_space.h"

#include "linear_solvers/reorderer.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "linear_solvers/scaling_solver.h"
#include "linear_solvers/mixedup_linear_solver.h"

#include "linear_solvers/preconditioner.h"
#include "linear_solvers/diagonal_preconditioner.h"
#include "linear_solvers/ilu0_preconditioner.h"
#include "linear_solvers/ilu_preconditioner.h"
//#include "linear_solvers/superlu_solver.h"
#include "linear_solvers/power_iteration_eigenvalue_solver.h"
#include "linear_solvers/deflated_gmres_solver.h"



namespace Kratos
{

namespace Python
{

template<typename TSpaceType, typename TLocalSpaceType>
void AddLinearSolversToPython_(const std::string& Prefix)
{
    typedef TSpaceType SpaceType;
    typedef typename SpaceType::DataType DataType;
    typedef typename SpaceType::ValueType ValueType;
    typedef TLocalSpaceType LocalSpaceType;

    typedef LinearSolver<SpaceType, LocalSpaceType> LinearSolverType;
    typedef IterativeSolver<SpaceType, LocalSpaceType> IterativeSolverType;
    typedef CGSolver<SpaceType, LocalSpaceType> CGSolverType;
    typedef DeflatedCGSolver<SpaceType, LocalSpaceType> DeflatedCGSolverType;
    typedef MixedUPLinearSolver<SpaceType, LocalSpaceType> MixedUPLinearSolverType;
    typedef BICGSTABSolver<SpaceType, LocalSpaceType> BICGSTABSolverType;
    typedef TFQMRSolver<SpaceType, LocalSpaceType> TFQMRSolverType;
    typedef ScalingSolver<SpaceType, LocalSpaceType> ScalingSolverType;
    typedef PowerIterationEigenvalueSolver<SpaceType, LocalSpaceType, LinearSolverType> PowerIterationEigenvalueSolverType;
    typedef DeflatedGMRESSolver<SpaceType, LocalSpaceType> DeflatedGMRESSolverType;

    bool (LinearSolverType::*pointer_to_solve)(typename LinearSolverType::SparseMatrixType& rA,
            typename LinearSolverType::VectorType& rX, typename LinearSolverType::VectorType& rB) = &LinearSolverType::Solve;

    using namespace boost::python;

    //****************************************************************************************************
    //preconditioners
    //****************************************************************************************************
    typedef Preconditioner<SpaceType, LocalSpaceType> PreconditionerType;

    class_<PreconditionerType, typename PreconditionerType::Pointer>((Prefix + "Preconditioner").c_str())
    .def(self_ns::str(self))
    ;

    typedef DiagonalPreconditioner<SpaceType, LocalSpaceType> DiagonalPreconditionerType;
    class_<DiagonalPreconditionerType, typename DiagonalPreconditionerType::Pointer, bases<PreconditionerType> >((Prefix + "DiagonalPreconditioner").c_str())
    .def(self_ns::str(self))
    ;

    typedef ILUPreconditioner<SpaceType, LocalSpaceType> ILUPreconditionerType;
    class_<ILUPreconditionerType, typename ILUPreconditionerType::Pointer, bases<PreconditionerType> >((Prefix + "ILUPreconditioner").c_str())
    .def(self_ns::str(self))
    ;

    typedef ILU0Preconditioner<SpaceType, LocalSpaceType> ILU0PreconditionerType;
    class_<ILU0PreconditionerType, typename ILU0PreconditionerType::Pointer, bases<PreconditionerType> >((Prefix + "ILU0Preconditioner").c_str())
    .def(self_ns::str(self))
    ;

    //****************************************************************************************************
    //linear solvers
    //****************************************************************************************************
    class_<LinearSolverType, typename LinearSolverType::Pointer>((Prefix + "LinearSolver").c_str())
    .def("Initialize",&LinearSolverType::Initialize)
    .def("Solve",pointer_to_solve)
    .def("Clear",&LinearSolverType::Clear)
    .def("SetEchoLevel", &LinearSolverType::SetEchoLevel)
    //.def("",&LinearSolverType::)
    .def(self_ns::str(self))
    ;

    class_<IterativeSolverType, typename IterativeSolverType::Pointer, bases<LinearSolverType> >((Prefix + "IterativeSolver").c_str())
    .def(self_ns::str(self))
    ;

    class_<CGSolverType, typename CGSolverType::Pointer, bases<IterativeSolverType> >((Prefix + "CGSolver").c_str())
    .def(init<ValueType>())
    .def(init<ValueType, unsigned int>())
    .def(init<ValueType, unsigned int, typename PreconditionerType::Pointer>())
    //.def("",&LinearSolverType::)
    .def(self_ns::str(self))
    ;

    class_<BICGSTABSolverType, typename BICGSTABSolverType::Pointer, bases<IterativeSolverType> >((Prefix + "BICGSTABSolver").c_str())
    .def(init<ValueType>())
    .def(init<ValueType, unsigned int>())
    .def(self_ns::str(self))
    .def(init<ValueType, unsigned int, typename PreconditionerType::Pointer>())
    .def("SetTolerance",&BICGSTABSolverType::SetTolerance)
    ;

    class_<TFQMRSolverType, typename TFQMRSolverType::Pointer, bases<IterativeSolverType> >((Prefix + "TFQMRSolver").c_str())
    .def(init<ValueType>())
    .def(init<ValueType, unsigned int>())
    .def(self_ns::str(self))
    .def(init<ValueType, unsigned int, typename PreconditionerType::Pointer>())
    ;

    class_<ScalingSolverType, typename ScalingSolverType::Pointer, bases<LinearSolverType> >((Prefix + "ScalingSolver").c_str())
    .def(init<typename LinearSolverType::Pointer, bool >())
    ;

    class_<PowerIterationEigenvalueSolverType, typename PowerIterationEigenvalueSolverType::Pointer, bases<LinearSolverType> >((Prefix + "PowerIterationEigenvalueSolver").c_str())
    .def(init<ValueType, unsigned int, unsigned int, typename LinearSolverType::Pointer>())
    ;

    typedef Reorderer<SpaceType, LocalSpaceType > ReordererType;
    typedef DirectSolver<SpaceType, LocalSpaceType, ReordererType > DirectSolverType;
    typedef SkylineLUFactorizationSolver<SpaceType, LocalSpaceType, ReordererType > SkylineLUFactorizationSolverType;

    class_<ReordererType, typename ReordererType::Pointer >((Prefix + "Reorderer").c_str())
    .def( init< >() )
    .def(self_ns::str(self))
    .def( "Initialize",&ReordererType::Initialize)
    .def( "Reorder",&ReordererType::Reorder)
    .def( "InverseReorder",&ReordererType::InverseReorder)
    ;

    class_<DirectSolverType, typename DirectSolverType::Pointer, bases<LinearSolverType> >((Prefix + "DirectSolver").c_str())
    .def( init< >() )
    .def(self_ns::str(self))
    ;

    class_<SkylineLUFactorizationSolverType, typename SkylineLUFactorizationSolverType::Pointer, bases<DirectSolverType> >((Prefix + "SkylineLUFactorizationSolver").c_str())
    .def(init< >())
    .def(self_ns::str(self))
    ;

    class_<DeflatedCGSolverType, typename DeflatedCGSolverType::Pointer, bases<IterativeSolverType> >((Prefix + "DeflatedCGSolver").c_str())
    .def(init<ValueType,bool,int>())
    .def(init<ValueType, unsigned int,bool,int>())
    .def(init<ValueType, unsigned int, typename PreconditionerType::Pointer,bool,int>())
    // .def(init<ValueType, unsigned int,  PreconditionerType::Pointer, ModelPart::Pointer>())
    // .def("",&LinearSolverType::)
    .def(self_ns::str(self))
    ;

    class_<MixedUPLinearSolverType, typename MixedUPLinearSolverType::Pointer, bases<IterativeSolverType> >((Prefix + "MixedUPLinearSolver").c_str(), init<typename LinearSolverType::Pointer, typename LinearSolverType::Pointer, ValueType, unsigned int, unsigned int >())
    .def(self_ns::str(self))
    ;

    class_<DeflatedGMRESSolverType, typename DeflatedGMRESSolverType::Pointer, bases<IterativeSolverType> >((Prefix + "DeflatedGMRESSolver").c_str(), init<typename LinearSolverType::Pointer, ValueType, unsigned int, unsigned int, unsigned int >())
    .def(self_ns::str(self))
    ;
}

void AddLinearSolversToPython()
{
    typedef KRATOS_DOUBLE_TYPE DataType;
    typedef UblasSpace<DataType, CompressedMatrix, Vector> SpaceType;
    typedef UblasSpace<DataType, Matrix, Vector> LocalSpaceType;

    AddLinearSolversToPython_<SpaceType, LocalSpaceType>("");

    typedef KRATOS_COMPLEX_TYPE ComplexType;
    typedef UblasSpace<ComplexType, ComplexCompressedMatrix, ComplexVector> ComplexSpaceType;
    typedef UblasSpace<ComplexType, ComplexMatrix, ComplexVector> LocalComplexSpaceType;

    AddLinearSolversToPython_<ComplexSpaceType, LocalComplexSpaceType>("Complex");
}

}  // namespace Python.

} // Namespace Kratos
