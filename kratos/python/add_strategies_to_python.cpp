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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


// Project includes
#include "includes/define.h"
#include "python/add_strategies_to_python.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy.h"
#include "solving_strategies/strategies/explicit_strategy.h"
//#include "solving_strategies/strategies/residualbased_arc_lenght_strategy.h"

//schemes
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/incremental_displacement_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "solving_strategies/convergencecriterias/or_criteria.h"
//#include "solving_strategies/convergencecriterias/and_criteria.h"

//Builder And Solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_deactivation.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_constraints.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_constraints_elementwise.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_constraints_deactivation.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_constraints_deactivation_elementwise.h"


//linear solvers
#include "linear_solvers/linear_solver.h"

//utilities
#include "python/pointer_vector_set_python_interface.h"


namespace Kratos
{

    namespace Python
    {
        using namespace boost::python;

        typedef KRATOS_DOUBLE_TYPE DataType;
        typedef UblasSpace<DataType, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<DataType, Matrix, Vector> LocalSpaceType;

        //ADDED BY PAOLO (next two)

        DataType Dot(SparseSpaceType& dummy, const SparseSpaceType::VectorType& rX, const SparseSpaceType::VectorType& rY)
        {
            return dummy.Dot(rX, rY);
        }

        void ScaleAndAdd(SparseSpaceType& dummy, const DataType A, const SparseSpaceType::VectorType& rX, const DataType B, SparseSpaceType::VectorType& rY)
        //(const DataType A,const  VectorType& rX, const DataType B, VectorType& rY) // rY = (A * rX) + (B * rY)
        {
            dummy.ScaleAndAdd(A, rX, B, rY);
        }

        void Mult(SparseSpaceType& dummy, SparseSpaceType::MatrixType& rA, SparseSpaceType::VectorType& rX, SparseSpaceType::VectorType& rY)
        //rY=A*rX (the product is stored inside the rY)
        {
            dummy.Mult(rA, rX, rY);
        }

        void TransposeMult(SparseSpaceType& dummy, SparseSpaceType::MatrixType& rA, SparseSpaceType::VectorType& rX, SparseSpaceType::VectorType& rY)
        //rY=A*rX (the product is stored inside the rY)
        {
            dummy.TransposeMult(rA, rX, rY);
        }

        SparseSpaceType::IndexType Size(SparseSpaceType& dummy, SparseSpaceType::VectorType const& rV)
        {
            return rV.size();
        }

        SparseSpaceType::IndexType Size1(SparseSpaceType& dummy, SparseSpaceType::MatrixType const& rM)
        {
            return rM.size1();
        }

        SparseSpaceType::IndexType Size2(SparseSpaceType& dummy, SparseSpaceType::MatrixType const& rM)
        {
            return rM.size2();
        }

        void ResizeMatrix(SparseSpaceType& dummy, SparseSpaceType::MatrixType& A, unsigned int i1, unsigned int i2)
        {
            dummy.Resize(A, i1, i2);
        }

        void ResizeVector(SparseSpaceType& dummy, SparseSpaceType::VectorType& x, unsigned int i1)
        {
            dummy.Resize(x, i1);
        }

        void SetToZeroMatrix(SparseSpaceType& dummy, SparseSpaceType::MatrixType& A)
        {
            dummy.SetToZero(A);
        }

        void SetToZeroVector(SparseSpaceType& dummy, SparseSpaceType::VectorType& x)
        {
            dummy.SetToZero(x);
        }

        void ClearMatrix(SparseSpaceType& dummy, SparseSpaceType::MatrixPointerType& A)
        {
            dummy.Clear(A);
        }

        void ClearVector(SparseSpaceType& dummy, SparseSpaceType::VectorPointerType& x)
        {
            dummy.Clear(x);
        }

        DataType TwoNorm(SparseSpaceType& dummy, const SparseSpaceType::VectorType& x)
        {
            return dummy.TwoNorm(x);
        }

        void UnaliasedAdd(SparseSpaceType& dummy, SparseSpaceType::VectorType& x, const DataType A, const SparseSpaceType::VectorType& rY) // x+= a*Y
        {
            dummy.UnaliasedAdd(x, A, rY);
        }

        void CopyVector(SparseSpaceType& dummy, SparseSpaceType::VectorType const& rX, SparseSpaceType::VectorType& rY)
        {
            dummy.Copy(rX, rY);
        }

        void CopyMatrix(SparseSpaceType& dummy, SparseSpaceType::MatrixType const& rX, SparseSpaceType::MatrixType& rY)
        {
            dummy.Copy(rX, rY);
        }

        void MoveMesh(Scheme< SparseSpaceType, LocalSpaceType >& dummy, ModelPart::NodesContainerType& rNodes)
        {
            for (ModelPart::NodeIterator i = rNodes.begin(); i != rNodes.end(); ++i)
            {
                const array_1d<DataType, 3> & disp = i->FastGetSolutionStepValue(DISPLACEMENT);
                (i)->X() = (i)->X0() + disp[0];
                (i)->Y() = (i)->Y0() + disp[1];
                (i)->Z() = (i)->Z0() + disp[2];
            }
        }

        SparseSpaceType::MatrixPointerType CreateEmptyMatrixPointer(SparseSpaceType& dummy)
        {
            return dummy.CreateEmptyMatrixPointer();
        }

        SparseSpaceType::VectorPointerType CreateEmptyVectorPointer(SparseSpaceType& dummy)
        {
            return dummy.CreateEmptyVectorPointer();
        }

        void PrintMatrixInfo(SparseSpaceType& dummy, const typename SparseSpaceType::MatrixType& rA, const int level)
        {
            dummy.PrintMatrixInfo(std::cout, rA, level);
        }

        void PrintVectorInfo(SparseSpaceType& dummy, const typename SparseSpaceType::VectorType& rX, const int level)
        {
            dummy.PrintVectorInfo(std::cout, rX, level);
        }

        //  boost::shared_ptr< CompressedMatrix > CreateEmptyMatrixPointer()
        //  {
        //      boost::shared_ptr<CompressedMatrix> pNewMat = boost::shared_ptr<CompressedMatrix>(new CompressedMatrix() );
        //      return pNewMat;
        //  }
        //
        //  boost::shared_ptr< Vector > CreateEmptyVectorPointer()
        //  {
        //      boost::shared_ptr<Vector > pNewVec = boost::shared_ptr<Vector >(new Vector() );
        //      return pNewVec;
        //  }

        CompressedMatrix& GetMatRef(boost::shared_ptr<CompressedMatrix>& dummy)
        {
            return *dummy;
        }

        Vector& GetVecRef(boost::shared_ptr<Vector>& dummy)
        {
            return *dummy;
        }

        template< class TOtherMatrixType >
        bool WriteMatrixMarketMatrix(SparseSpaceType& dummy, const char *FileName, TOtherMatrixType &M, bool Symmetric)
        {
            return dummy.WriteMatrixMarketMatrix(FileName, M, Symmetric);
        }

        template< class VectorType >
        bool WriteMatrixMarketVector(SparseSpaceType& dummy, const char *FileName, VectorType& V)
        {
            return dummy.WriteMatrixMarketVector(FileName, V);
        }

        template< class CompressedMatrixType >
        bool ReadMatrixMarketMatrix_(SparseSpaceType& dummy, const char *FileName, CompressedMatrixType& M)
        {
            return ReadMatrixMarketMatrix<CompressedMatrixType>(FileName, M);
        }

        template< class VectorType >
        bool ReadMatrixMarketVector_(SparseSpaceType& dummy, const char *FileName, VectorType& M)
        {
            return ReadMatrixMarketVector<VectorType>(FileName, M);
        }

        void AddStrategiesToPython()
        {
            class_< boost::shared_ptr<CompressedMatrix> >("CompressedMatrixPointer", init<boost::shared_ptr<CompressedMatrix> >())
                    .def("GetReference", GetMatRef, return_value_policy<reference_existing_object > ())
                    //                  .def("GetReference", GetRef, return_internal_reference<1>() )
                    ;

            class_< boost::shared_ptr<Vector> >("VectorPointer", init< boost::shared_ptr<Vector> >())
                    .def("GetReference", GetVecRef, return_value_policy<reference_existing_object > ())
                    ;

            typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
            typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
            typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
            typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ::Pointer TConvergenceCriteriaPointer;

            //********************************************************************
            //********************************************************************
            //strategy base class
            class_< BaseSolvingStrategyType, boost::noncopyable > ("SolvingStrategy", init < ModelPart&, bool >())
                    .def("Predict", &BaseSolvingStrategyType::Predict)
                    .def("Initialize", &BaseSolvingStrategyType::Initialize)
                    .def("Solve", &BaseSolvingStrategyType::Solve)
                    .def("IsConverged", &BaseSolvingStrategyType::IsConverged)
                    .def("CalculateOutputData", &BaseSolvingStrategyType::CalculateOutputData)
                    .def("SetEchoLevel", &BaseSolvingStrategyType::SetEchoLevel)
                    .def("GetEchoLevel", &BaseSolvingStrategyType::GetEchoLevel)
                    .def("SetRebuildLevel", &BaseSolvingStrategyType::SetRebuildLevel)
                    .def("GetRebuildLevel", &BaseSolvingStrategyType::GetRebuildLevel)
                    .def("SetMoveMeshFlag", &BaseSolvingStrategyType::SetMoveMeshFlag)
                    .def("MoveMeshFlag", &BaseSolvingStrategyType::MoveMeshFlag)
                    .def("MoveMesh", &BaseSolvingStrategyType::MoveMesh)
                    .def("Clear", &BaseSolvingStrategyType::Clear)
                    .def("Check", &BaseSolvingStrategyType::Check)
                    .def("InitializeSolutionStep", &BaseSolvingStrategyType::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &BaseSolvingStrategyType::FinalizeSolutionStep)
                    .def("SolveSolutionStep", &BaseSolvingStrategyType::SolveSolutionStep)
                    //.def("GetModelPart", &BaseSolvingStrategyType::GetModelPart )
                    ;

            typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > TConvergenceCriteriaType;
            typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;

            class_< ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, bases< BaseSolvingStrategyType >, boost::noncopyable >
                    ("ResidualBasedLinearStrategy",init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, bool, bool, bool, bool >())
                    .def(init < ModelPart& ,  BaseSchemeType::Pointer, LinearSolverType::Pointer, BuilderAndSolverType::Pointer, bool, bool, bool,  bool  >())
                    .def("GetResidualNorm", &ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetResidualNorm)
                    .def("SetBuilderAndSolver", &ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetBuilderAndSolver)
                    ;


            class_< ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, bases< BaseSolvingStrategyType >, boost::noncopyable >
                    ("ResidualBasedNewtonRaphsonStrategy", init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, bool, bool, bool >())
                    .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >())
                    .def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetMaxIterationNumber)
                    .def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetMaxIterationNumber)
                    .def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetKeepSystemConstantDuringIterations)
                    .def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetKeepSystemConstantDuringIterations)
                    ;

            class_< AdaptiveResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >, bases< BaseSolvingStrategyType >, boost::noncopyable >
                    ("AdaptiveResidualBasedNewtonRaphsonStrategy",
                    init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, TConvergenceCriteriaType::Pointer, int, int, bool, bool, bool, DataType, DataType, int
                    >())
                    ;

            class_< ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
                    bases< BaseSolvingStrategyType >,  boost::noncopyable >
                    ("Explicit_Strategy",
                    init<ModelPart&, int, bool >() )
                    //AssembleLoop loops the elements calling AddExplicitContribution. Using processinfo the element is the one who "decides" which variable to modify.
                    .def("AssembleLoop",&ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::AssembleLoop)
                    //once the assembleloop has been performed, the variable must be normalized. (for example with the nodal mass or the nodal area). Loop on nodes.
                    .def("NormalizeVariable",&ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::NormalizeVariable)
                    //ExplicitUpdateLoop modifies a vectorial variable by adding another variable (the RHS, PRESS_PROJ,etc) multiplied by a user-given factor (ie delta_time)
                    .def("ExplicitUpdateLoop",&ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::ExplicitUpdateLoop)
                    //initialize and finalize.
                    .def("InitializeSolutionStep",&ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::InitializeSolutionStep)
                    .def("FinalizeSolutionStep",&ExplicitStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::FinalizeSolutionStep)
                    ;

            //********************************************************************
            //********************************************************************
            class_< BaseSchemeType, boost::noncopyable >
                    ("Scheme", init< >())
                    .def("Initialize", &BaseSchemeType::Initialize)
                    .def("SchemeIsInitialized", &BaseSchemeType::SchemeIsInitialized)
                    .def("ElementsAreInitialized", &BaseSchemeType::ElementsAreInitialized)
                    .def("ConditionsAreInitialized", &BaseSchemeType::ConditionsAreInitialized)
                    .def("InitializeElements", &BaseSchemeType::InitializeElements)
                    .def("InitializeConditions", &BaseSchemeType::InitializeConditions)
                    .def("InitializeSolutionStep", &BaseSchemeType::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &BaseSchemeType::FinalizeSolutionStep)
                    .def("InitializeNonLinIteration", &BaseSchemeType::InitializeNonLinIteration)
                    .def("FinalizeNonLinIteration", &BaseSchemeType::FinalizeNonLinIteration)
                    .def("Predict", &BaseSchemeType::Predict)
                    .def("Update", &BaseSchemeType::Update)
                    .def("CalculateOutputData", &BaseSchemeType::CalculateOutputData)
                    .def("Clean", &BaseSchemeType::Clean)
                    .def("Clear",&BaseSchemeType::Clear)
                    .def("MoveMesh", MoveMesh)
                    .def("Check", &BaseSchemeType::Check)
                    ;

            class_< ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType>,
                    bases< BaseSchemeType >, boost::noncopyable >
                    (
                    "ResidualBasedIncrementalUpdateStaticScheme", init< >()
                    );

            class_< ResidualBasedIncrementalUpdateStaticSchemeSlip< SparseSpaceType, LocalSpaceType>,
                    bases< ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType> >,
                    boost::noncopyable >
                    ("ResidualBasedIncrementalUpdateStaticSchemeSlip", init<unsigned int, unsigned int>());


            //********************************************************************
            //********************************************************************
            //convergence criteria base class
            class_< ConvergenceCriteria< SparseSpaceType, LocalSpaceType >, boost::noncopyable > ("ConvergenceCriteria", init<>())
                    .def("SetActualizeRHSFlag", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::SetActualizeRHSFlag)
                    .def("GetActualizeRHSflag", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::GetActualizeRHSflag)
                    .def("PreCriteria", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::PreCriteria)
                    .def("PostCriteria", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::PostCriteria)
                    .def("Initialize", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::Initialize)
                    .def("InitializeNonLinearIteration", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::InitializeNonLinearIteration)
                    .def("InitializeSolutionStep", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::InitializeSolutionStep)
                    .def("FinalizeNonLinearIteration", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::FinalizeNonLinearIteration)
                    .def("FinalizeSolutionStep", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::FinalizeSolutionStep)
                    .def("Check", &ConvergenceCriteria<SparseSpaceType, LocalSpaceType >::Check)
                    ;

            class_< DisplacementCriteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("DisplacementCriteria", init< DataType, DataType>())
                    .def("SetEchoLevel", &DisplacementCriteria<SparseSpaceType, LocalSpaceType >::SetEchoLevel);

            class_< IncrementalDisplacementCriteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("IncrementalDisplacementCriteria", init< DataType, DataType>());

            class_<ResidualCriteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("ResidualCriteria", init< DataType, DataType>());

            /*          class_< ResidualCriteria< SparseSpaceType >,
                                             bases<ConvergenceCriteria< SparseSpaceType > >,
                                             boost::noncopyable >
                                            ("ResidualCriteria", init<Model::Pointer, DataType >() );

                                    class_< AndCriteria< SparseSpaceType >,
                                             bases<ConvergenceCriteria< SparseSpaceType > >,
                                             boost::noncopyable >
                                            ("AndCriteria", init<Model::Pointer, ConvergenceCriteria< SparseSpaceType >::Pointer, ConvergenceCriteria< SparseSpaceType >::Pointer >()*/
            //);


            class_<And_Criteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("AndCriteria", init<TConvergenceCriteriaPointer, TConvergenceCriteriaPointer > ());

            class_<Or_Criteria<SparseSpaceType, LocalSpaceType >,
                    bases<ConvergenceCriteria< SparseSpaceType, LocalSpaceType > >,
                    boost::noncopyable >
                    ("OrCriteria", init<TConvergenceCriteriaPointer, TConvergenceCriteriaPointer > ());

            //********************************************************************
            //********************************************************************
            //

            //Builder and Solver
            class_< BuilderAndSolverType::DofsArrayType, boost::noncopyable > ("DofsArrayType", init<>());
            BuilderAndSolverType::DofsArrayType&(BuilderAndSolverType::*BuilderAndSolver_GetDofSet)() = &BuilderAndSolverType::GetDofSet;

            class_< BuilderAndSolverType, boost::noncopyable > ("BuilderAndSolver", init<LinearSolverType::Pointer > ())
                    .def("SetCalculateReactionsFlag", &BuilderAndSolverType::SetCalculateReactionsFlag)
                    .def("GetCalculateReactionsFlag", &BuilderAndSolverType::GetCalculateReactionsFlag)
                    .def("SetDofSetIsInitializedFlag", &BuilderAndSolverType::SetDofSetIsInitializedFlag)
                    .def("GetDofSetIsInitializedFlag", &BuilderAndSolverType::GetDofSetIsInitializedFlag)
                    .def("SetReshapeMatrixFlag", &BuilderAndSolverType::SetReshapeMatrixFlag)
                    .def("GetReshapeMatrixFlag", &BuilderAndSolverType::GetReshapeMatrixFlag)
                    .def("GetEquationSystemSize", &BuilderAndSolverType::GetEquationSystemSize)
                    .def("BuildLHS", &BuilderAndSolverType::BuildLHS)
                    .def("BuildRHS", &BuilderAndSolverType::BuildRHS)
                    .def("Build", &BuilderAndSolverType::Build)
                    .def("SystemSolve", &BuilderAndSolverType::SystemSolve)
                    .def("BuildAndSolve", &BuilderAndSolverType::BuildAndSolve)
                    .def("BuildRHSAndSolve", &BuilderAndSolverType::BuildRHSAndSolve)
                    .def("ApplyDirichletConditions", &BuilderAndSolverType::ApplyDirichletConditions)
                    .def("SetUpDofSet", &BuilderAndSolverType::SetUpDofSet)
                    .def("GetDofSet", BuilderAndSolver_GetDofSet, return_internal_reference<>())
                    .def("SetUpSystem", &BuilderAndSolverType::SetUpSystem)
                    .def("ResizeAndInitializeVectors", &BuilderAndSolver_ResizeAndInitializeVectors1<BuilderAndSolverType>)
                    .def("ResizeAndInitializeVectors", &BuilderAndSolver_ResizeAndInitializeVectors2<BuilderAndSolverType>)
                    .def("InitializeSolutionStep", &BuilderAndSolverType::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &BuilderAndSolverType::FinalizeSolutionStep)
                    .def("CalculateReactions", &BuilderAndSolverType::CalculateReactions)
                    .def("Clear", &BuilderAndSolverType::Clear)
                    .def("Check", &BuilderAndSolverType::Check)
                    .def("SetEchoLevel", &BuilderAndSolverType::SetEchoLevel)
                    .def("GetEchoLevel", &BuilderAndSolverType::GetEchoLevel)
                    ;

            typedef ResidualBasedEliminationBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverType;
            class_< ResidualBasedEliminationBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > ("ResidualBasedEliminationBuilderAndSolver", init< LinearSolverType::Pointer > ());

            typedef ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverType;
            class_< ResidualBasedBlockBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > ("ResidualBasedBlockBuilderAndSolver", init< LinearSolverType::Pointer > ());


            typedef ResidualBasedEliminationBuilderAndSolverDeactivation< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedEliminationBuilderAndSolverDeactivationType;
            class_< ResidualBasedEliminationBuilderAndSolverDeactivationType, bases<BuilderAndSolverType>, boost::noncopyable > ("ResidualBasedEliminationBuilderAndSolverDeactivation", init< LinearSolverType::Pointer > ());

            typedef ResidualBasedBlockBuilderAndSolverWithConstraints< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverWithConstraintsType;
            class_< ResidualBasedBlockBuilderAndSolverWithConstraintsType, bases<BuilderAndSolverType>, boost::noncopyable > ("ResidualBasedBlockBuilderAndSolverWithConstraints", init< LinearSolverType::Pointer > ());

            typedef ResidualBasedBlockBuilderAndSolverWithConstraintsElementWise< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverWithConstraintsElementWiseType;
            class_< ResidualBasedBlockBuilderAndSolverWithConstraintsElementWiseType, bases<ResidualBasedBlockBuilderAndSolverWithConstraintsType>, boost::noncopyable > ("ResidualBasedBlockBuilderAndSolverWithConstraintsElementWise", init< LinearSolverType::Pointer > ());

            typedef ResidualBasedBlockBuilderAndSolverWithConstraintsDeactivation< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverWithConstraintsDeactivationType;
            class_< ResidualBasedBlockBuilderAndSolverWithConstraintsDeactivationType, bases<BuilderAndSolverType>, boost::noncopyable > ("ResidualBasedBlockBuilderAndSolverWithConstraintsDeactivation", init< LinearSolverType::Pointer > ());

            typedef ResidualBasedBlockBuilderAndSolverWithConstraintsDeactivationElementWise< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverWithConstraintsDeactivationElementWiseType;
            class_< ResidualBasedBlockBuilderAndSolverWithConstraintsDeactivationElementWiseType, bases<BuilderAndSolverType>, boost::noncopyable > ("ResidualBasedBlockBuilderAndSolverWithConstraintsDeactivationElementWise", init< LinearSolverType::Pointer > ());

            //
            // class_< ResidualBasedBlockBuilderAndSolverWithConstraintsType,
            //             ResidualBasedBlockBuilderAndSolverWithConstraintsType::Pointer,
            //             BuilderAndSolverType>
            //             (m,"ResidualBasedBlockBuilderAndSolverWithConstraints").def(init< LinearSolverType::Pointer > ());
            //********************************************************************
            //********************************************************************

            class_< SparseSpaceType, boost::noncopyable > ("UblasSparseSpace", init<>())
                    .def("ClearMatrix", ClearMatrix)
                    .def("ClearVector", ClearVector)
                    .def("ResizeMatrix", ResizeMatrix)
                    .def("ResizeVector", ResizeVector)
                    .def("SetToZeroMatrix", SetToZeroMatrix)
                    .def("SetToZeroVector", SetToZeroVector)
                    .def("TwoNorm", TwoNorm)
                    .def("Copy", CopyVector)
                    .def("Copy", CopyMatrix)
                    //the dot product of two vectors
                    .def("Dot", Dot)
                    //the matrix-vector multiplication
                    .def("Mult", Mult)
                    .def("TransposeMult", TransposeMult)
                    .def("Size", Size)
                    .def("Size1", Size1)
                    .def("Size2", Size2)
                    .def("UnaliasedAdd", UnaliasedAdd)
                    .def("ScaleAndAdd", ScaleAndAdd)
                    .def("CreateEmptyMatrixPointer", CreateEmptyMatrixPointer)
                    .def("CreateEmptyVectorPointer", CreateEmptyVectorPointer)
                    .def("WriteMatrixMarketMatrix", WriteMatrixMarketMatrix<CompressedMatrix>)
                    .def("WriteMatrixMarketVector", WriteMatrixMarketVector<Vector>)
                    .def("ReadMatrixMarketMatrix", ReadMatrixMarketMatrix_<CompressedMatrix>)
                    .def("ReadMatrixMarketVector", ReadMatrixMarketVector_<Vector>)
                    .def("PrintMatrixInfo", PrintMatrixInfo)
                    .def("PrintVectorInfo", PrintVectorInfo)
                    .def(self_ns::str(self))
                    ;
        }

    } // namespace Python.

} // Namespace Kratos
