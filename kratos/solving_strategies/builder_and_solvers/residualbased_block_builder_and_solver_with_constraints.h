//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborators:   Vicente Mataix
//
//
#if !defined(KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS )
#define  KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS


/* System includes */
#include <iostream>
#include <unordered_set>
#include <unordered_map>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/model_part.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"
#include "includes/kratos_flags.h"
#include "utilities/sparse_matrix_multiplication_utility.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class ResidualBasedEliminationBuilderAndSolver
 * @ingroup KratosCore
 * @brief Current class provides an implementation for standard builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual)
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information.
 * Calculation of the reactions involves a cost very similiar to the calculation of the total residual
 * @author Riccardo Rossi
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver, //= LinearSolver<TSparseSpace,TDenseSpace>
         class TModelPartType
         >
class ResidualBasedBlockBuilderAndSolverWithConstraints
    : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver, TModelPartType >
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverWithConstraints);

    /// Definition of the base class
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver, TModelPartType> BaseType;

    /// Definition of the classes from the base class
    typedef typename BaseType::ModelPartType ModelPartType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::ValueType ValueType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::NodesContainerType NodesContainerType;
    typedef typename BaseType::ElementsContainerType ElementsContainerType;
    typedef typename BaseType::ConditionsContainerType ConditionsContainerType;

    /// Additional definitions
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::NodeType NodeType;
    typedef typename BaseType::ElementType ElementType;
    typedef typename BaseType::ConditionType ConditionType;
    typedef typename ElementType::EquationIdVectorType EquationIdVectorType;
    typedef typename ElementType::DofsVectorType DofsVectorType;
    typedef typename MatrixVectorTypeSelector<TDataType>::CompressedMatrixType CompressedMatrixType;

    /// DoF types definition
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (with parameters)
     */
    // explicit ResidualBasedBlockBuilderAndSolverWithConstraints(
    //     typename TLinearSolver::Pointer pNewLinearSystemSolver,
    //     Parameters ThisParameters
    //     ) : BaseType(pNewLinearSystemSolver)
    // {
    //     // Validate default parameters
    //     Parameters default_parameters = Parameters(R"(
    //     {
    //         "name" : "ResidualBasedBlockBuilderAndSolverWithConstraints"
    //     })" );
    //
    //     ThisParameters.ValidateAndAssignDefaults(default_parameters);
    // }

    /**
     * @brief Default constructor.
     */
    explicit ResidualBasedBlockBuilderAndSolverWithConstraints(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    ~ResidualBasedBlockBuilderAndSolverWithConstraints() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to perform the build of the RHS. The vector could be sized as the total number
     * of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param b The RHS vector
     */
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPartType& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& b) override
    {
        KRATOS_TRY
        if(!pScheme)
        {
            KRATOS_ERROR << "No scheme provided!";
        }
        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        std::cout << "Total number of elements in assembly: " << nelements << std::endl;
        std::cout << "Total number of conditions in assembly: " << nconditions << std::endl;
        std::cout << "Total number of constraints in assembly: " << rModelPart.NumberOfMasterSlaveConstraints() << std::endl;
        std::cout << "Number of threads: " << OpenMPUtils::GetNumThreads() << std::endl;

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        auto el_begin = rModelPart.ElementsBegin();
        auto cond_begin = rModelPart.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        typename ElementType::EquationIdVectorType EquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

        #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
        {
            # pragma omp for  schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++)
            {
                auto it = el_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    element_is_active = (it)->Is(ACTIVE);

                if (element_is_active)
                {
                    //calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);

                    // clean local elemental memory
                    pScheme->CleanMemory(*it);
                }
            }

            #pragma omp for  schedule(guided, 512)
            for (int k = 0; k < nconditions; k++)
            {
                auto it = cond_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool condition_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    condition_is_active = (it)->Is(ACTIVE);

                if (condition_is_active)
                {
                    //calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);

                    // clean local elemental memory
                    pScheme->CleanMemory(*it);
                }
            }
        }

        const double stop_build = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "Build time: " << stop_build - start_build << std::endl;
        }

        if (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "Finished parallel building" << std::endl;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building of the LHS
     * @details Depending on the implementation choosen the size of the matrix could
     * be equal to the total number of Dofs or to the number of unrestrained dofs
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     */
    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPartType& rModelPart,
        TSystemMatrixType& A) override
    {
        KRATOS_TRY

        TSystemVectorType tmp(A.size1(), 0.0);
        this->Build(pScheme, rModelPart, A, tmp);

        KRATOS_CATCH("")
    }

    /**
     * @brief Build a rectangular matrix of size n*N where "n" is the number of unrestrained degrees of freedom
     * and "N" is the total number of degrees of freedom involved.
     * @details This matrix is obtained by building the total matrix without the lines corresponding to the fixed
     * degrees of freedom (but keeping the columns!!)
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     */
    void BuildLHS_CompleteOnFreeRows(
        typename TSchemeType::Pointer pScheme,
        ModelPartType& rModelPart,
        TSystemMatrixType& A) override
    {
        KRATOS_TRY

        TSystemVectorType tmp(A.size1(), 0.0);
        this->Build(pScheme, rModelPart, A, tmp);

        KRATOS_CATCH("")
    }

    /**
     * @brief This is a call to the linear system solver
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void SystemSolve(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY

        double start_solve = OpenMPUtils::GetCurrentTime();

        ValueType norm_b;
        if (TSparseSpace::Size(b) != 0)
            norm_b = std::abs(TSparseSpace::TwoNorm(b));
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            //do solve
            BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
        }
        else
            TSparseSpace::SetToZero(Dx);

        if(mT.size1() != 0) //if there are master-slave constraints
        {
            //recover solution of the original problem
            TSystemVectorType Dxmodified = Dx;

            TSparseSpace::Mult(mT, Dxmodified, Dx);
        }

        double stop_solve = OpenMPUtils::GetCurrentTime();
        std::cout << "System Solve time: " << stop_solve - start_solve << "s" << std::endl;

        //prints informations about the current time
        if (this->GetEchoLevel() > 1)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << *(BaseType::mpLinearSystemSolver) << std::endl;
        }

        KRATOS_CATCH("")
    }

    void SystemSolveWithPhysics(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        ModelPartType& rModelPart
    )
    {
        if(rModelPart.MasterSlaveConstraints().size() != 0) {
            TSystemVectorType Dxmodified(b.size());

            InternalSystemSolveWithPhysics(A, Dxmodified, b, rModelPart);

            //recover solution of the original problem
            if (mMasterIds.size() != mT.size1())
            {
                double start_recover = OpenMPUtils::GetCurrentTime();
                TSparseSpace::Mult(mT, Dxmodified, Dx);
                double stop_recover = OpenMPUtils::GetCurrentTime();
                KRATOS_WATCH(norm_2(Dx))
                std::cout << "Recover back the system time: " << stop_recover - start_recover << "s" << std::endl;
            }
            else
            {
                TSparseSpace::Copy(Dxmodified, Dx);
                KRATOS_WATCH(norm_2(Dx))
            }
        } else {
            InternalSystemSolveWithPhysics(A, Dx, b, rModelPart);
        }
    }

    /**
      *@brief This is a call to the linear system solver (taking into account some physical particularities of the problem)
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    void InternalSystemSolveWithPhysics(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        ModelPartType& rModelPart
    )
    {
        KRATOS_TRY

        double start_solve = OpenMPUtils::GetCurrentTime();
        std::cout << "Begin Internal-System-Solve-With-Physics" << std::endl;

        ValueType norm_b;
        if (TSparseSpace::Size(b) != 0)
            norm_b = std::abs(TSparseSpace::TwoNorm(b));
        else
            norm_b = 0.00;

        if (norm_b != 0.00) {
            //provide physical data as needed
            if(BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded() )
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(A, Dx, b, BaseType::mDofSet, rModelPart);

            //do solve
            BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
        } else {
            TSparseSpace::SetToZero(Dx);
            std::cout << "ATTENTION! setting the RHS to zero!" << std::endl;
            // KRATOS_WARNING("ResidualBasedBlockBuilderAndSolverWithConstraints") << "ATTENTION! setting the RHS to zero!" << std::endl;
        }

        double stop_solve = OpenMPUtils::GetCurrentTime();
        std::cout << "Internal-System-Solve-With-Physics time: " << stop_solve - start_solve << "s" << std::endl;

        // Prints informations about the current time
        if (this->GetEchoLevel() > 1)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << *(BaseType::mpLinearSystemSolver) << std::endl;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building and solving phase at the same time.
     * @details It is ideally the fastest and safer function to use when it is possible to solve
     * just after building
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPartType& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        Timer::Start("Build");

        Build(pScheme, rModelPart, A, b);

        Timer::Stop("Build");

        if(rModelPart.MasterSlaveConstraints().size() != 0) {
            Timer::Start("ApplyConstraints");
            KRATOS_WATCH(mMasterIds.size())
            KRATOS_WATCH(mT.size1())
            KRATOS_WATCH(mT.size2())
            if (mMasterIds.size() != mT.size1())
                ApplyConstraints(pScheme,A,Dx,b,rModelPart);
            Timer::Stop("ApplyConstraints");
        }

        ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

        if ( this->GetEchoLevel() == 3)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "Before the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;
        }

        const double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");

        SystemSolveWithPhysics(A, Dx, b, rModelPart);

        Timer::Stop("Solve");
        const double stop_solve = OpenMPUtils::GetCurrentTime();

        if (this->GetEchoLevel() >=1 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "System solve time: " << stop_solve - start_solve << std::endl;
        }

        if ( this->GetEchoLevel() == 3)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "After the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Corresponds to the previews, but the System's matrix is considered already built and only the RHS is built again
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPartType& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        BuildRHS(pScheme, rModelPart, b);

        Timer::Stop("Build");

        if(rModelPart.MasterSlaveConstraints().size() != 0) {
            Timer::Start("ApplyRHSConstraints");
            KRATOS_WATCH(mMasterIds.size())
            KRATOS_WATCH(mT.size1())
            KRATOS_WATCH(mT.size2())
            if (mMasterIds.size() != mT.size1())
                ApplyRHSConstraints(pScheme, rModelPart, b);
            Timer::Stop("ApplyRHSConstraints");
        }

        ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

        if ( this->GetEchoLevel() == 3)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "Before the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;
        }

        const double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");

        SystemSolveWithPhysics(A, Dx, b, rModelPart);

        Timer::Stop("Solve");
        const double stop_solve = OpenMPUtils::GetCurrentTime();

        if (this->GetEchoLevel() >=1 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "System solve time: " << stop_solve - start_solve << std::endl;
        }

        if ( this->GetEchoLevel() == 3)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "After the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the build of the RHS.
     * @details The vector could be sized as the total number of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPartType& rModelPart,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        BuildRHSNoDirichlet(pScheme,rModelPart,b);

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(ndofs)
        for (int k = 0; k<ndofs; k++)
        {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;
            const std::size_t i = dof_iterator->EquationId();

            if (dof_iterator->IsFixed())
                b[i] = 0.0;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Builds the list of the DofSets involved in the problem by "asking" to each element
     * and condition its Dofs.
     * @details The list of dofs is stores insde the BuilderAndSolver as it is closely connected to the
     * way the matrix and RHS are built
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPartType& rModelPart
    ) override
    {
        KRATOS_TRY

        if ( this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "Setting up the dofs" << std::endl;
        }

        //Gets the array of elements from the modeler
        ElementsContainerType& r_elements_array = rModelPart.Elements();
        const int number_of_elements = static_cast<int>(r_elements_array.size());

        DofsVectorType dof_list, second_dof_list; // NOTE: The second dof list is only used on constraints to include master/slave relations

        unsigned int nthreads = OpenMPUtils::GetNumThreads();

        typedef std::set < typename DofType::Pointer >  set_type;

        if ( this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "Number of threads" << nthreads << "\n" << std::endl;
        }

        if ( this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "Initializing element loop" << std::endl;
        }

        /**
         * Here we declare three sets.
         * - The global set: Contains all the DoF of the system
         * - The slave set: The DoF that are not going to be solved, due to MPC formulation
         */
        set_type dof_global_set;

        #ifdef ENABLE_PARALLEL_SETUP_DOFSET
        #pragma omp parallel firstprivate(dof_list, second_dof_list)
        {
        #endif
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // We create the temporal set
            set_type dofs_tmp_set;

            // Gets the array of elements from the modeler
            #ifdef ENABLE_PARALLEL_SETUP_DOFSET
            #pragma omp for schedule(guided, 512) nowait
            #endif
            for (int i = 0; i < number_of_elements; ++i) {
                auto it_elem = r_elements_array.begin() + i;

                // Gets list of Dof involved on every element
                pScheme->GetDofList(*it_elem, dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
            }

            // Gets the array of conditions from the modeler
            ConditionsContainerType& r_conditions_array = rModelPart.Conditions();
            const int number_of_conditions = static_cast<int>(r_conditions_array.size());
            #ifdef ENABLE_PARALLEL_SETUP_DOFSET
            #pragma omp for  schedule(guided, 512) nowait
            #endif
            for (int i = 0; i < number_of_conditions; ++i) {
                auto it_cond = r_conditions_array.begin() + i;

                // Gets list of Dof involved on every element
                pScheme->GetDofList(*it_cond, dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
            }

            // Gets the array of constraints from the modeler
            auto& r_constraints_array = rModelPart.MasterSlaveConstraints();
            const int number_of_constraints = static_cast<int>(r_constraints_array.size());
            #ifdef ENABLE_PARALLEL_SETUP_DOFSET
            #pragma omp for  schedule(guided, 512) nowait
            #endif
            for (int i = 0; i < number_of_constraints; ++i) {
                auto it_const = r_constraints_array.begin() + i;

                // Gets list of Dof involved on every element
                it_const->GetDofList(dof_list, second_dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                dofs_tmp_set.insert(second_dof_list.begin(), second_dof_list.end());
            }

            // We merge all the sets in one thread
            #ifdef ENABLE_PARALLEL_SETUP_DOFSET
            #pragma omp critical
            {
            #endif
                dof_global_set.insert(dofs_tmp_set.begin(), dofs_tmp_set.end());
            #ifdef ENABLE_PARALLEL_SETUP_DOFSET
            }
            #endif
        #ifdef ENABLE_PARALLEL_SETUP_DOFSET
        }
        #endif

        if ( this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "Initializing ordered array filling\n" << std::endl;
        }

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        Doftemp.reserve(dof_global_set.size());
        for (auto it= dof_global_set.begin(); it!= dof_global_set.end(); it++)
        {
            Doftemp.push_back( *it );
        }
        Doftemp.Sort();

        BaseType::mDofSet = Doftemp;

        //Throws an exception if there are no Degrees Of Freedom involved in the analysis
        if(BaseType::mDofSet.size() == 0)
        {
            KRATOS_ERROR << "No degrees of freedom!";
        }

        if ( this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "Number of degrees of freedom:" << BaseType::mDofSet.size() << std::endl;
        }

        BaseType::mDofSetIsInitialized = true;

        if ( this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "Finished setting up the dofs" << std::endl;
        }

        if ( this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "End of setup dof set\n" << std::endl;
        }

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        // This is tobe done only in debug mode
        if (BaseType::GetCalculateReactionsFlag()) {
            for (auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator) {
                    if(!(dof_iterator->HasReaction()))
                    {
                        KRATOS_ERROR << "Reaction variable not set for the following Node: " << dof_iterator->Id() << ". Not possible to calculate reactions.";
                    }
            }
        }
#endif

        KRATOS_CATCH("");
    }

    /**
     * @brief Organises the dofset in order to speed up the building phase
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpSystem(
        ModelPartType& rModelPart
    ) override
    {
        //int free_id = 0;
        BaseType::mEquationSystemSize = BaseType::mDofSet.size();
        int ndofs = static_cast<int>(BaseType::mDofSet.size());

        #pragma omp parallel for firstprivate(ndofs)
        for (int i = 0; i < static_cast<int>(ndofs); i++) {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + i;
            dof_iterator->SetEquationId(i);

        }
    }

    //**************************************************************************
    //**************************************************************************

    void ResizeAndInitializeVectors(
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ModelPartType& rModelPart
    ) //override
    {
        KRATOS_TRY
        if (pA == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
            pA.swap(pNewA);
        }
        if (pDx == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
            pDx.swap(pNewDx);
        }
        if (pb == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
            pb.swap(pNewb);
        }
        if (BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(0) );
            BaseType::mpReactionsVector.swap(pNewReactionsVector);
        }

        TSystemMatrixType& A = *pA;
        TSystemVectorType& Dx = *pDx;
        TSystemVectorType& b = *pb;

        //resizing the system vectors and matrix
        if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        {
            A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
            ConstructMatrixStructure(A, rModelPart);
        }
        else
        {
            if (A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
            {
                KRATOS_ERROR << "The equation system size has changed during the simulation. This is not permited.";
                A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
                ConstructMatrixStructure(A, rModelPart);
            }
        }
        if (Dx.size() != BaseType::mEquationSystemSize)
            Dx.resize(BaseType::mEquationSystemSize, false);
        if (b.size() != BaseType::mEquationSystemSize)
            b.resize(BaseType::mEquationSystemSize, false);

        //if needed resize the vector for the calculation of reactions
        if(BaseType::mCalculateReactionsFlag == true)
        {
            unsigned int ReactionsVectorSize = BaseType::mDofSet.size()-BaseType::mEquationSystemSize;
            if(BaseType::mpReactionsVector->size() != ReactionsVectorSize)
                BaseType::mpReactionsVector->resize(ReactionsVectorSize,false);
        }

        ConstructMasterSlaveConstraintsStructure(rModelPart);

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    void InitializeSolutionStep(
        ModelPartType& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        KRATOS_TRY

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        // Getting process info
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // Computing constraints
        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        #pragma omp parallel for schedule(guided, 512) firstprivate(n_constraints, constraints_begin)
        for (int k = 0; k < n_constraints; ++k) {
            auto it = constraints_begin + k;
            it->InitializeSolutionStep(r_process_info); // Here each constraint constructs and stores its T and C matrices. Also its equation slave_ids.
        }

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    void FinalizeSolutionStep(
        ModelPartType& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        // Getting process info
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // Computing constraints
        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        const auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        #pragma omp parallel for schedule(guided, 512) firstprivate(n_constraints, constraints_begin)
        for (int k = 0; k < n_constraints; ++k) {
            auto it = constraints_begin + k;
            it->FinalizeSolutionStep(r_process_info);
        }
    }

    //**************************************************************************
    //**************************************************************************

    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPartType& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        TSparseSpace::SetToZero(b);

        //refresh RHS to have the correct reactions
        BuildRHSNoDirichlet(pScheme, rModelPart, b);

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(ndofs)
        for (int k = 0; k<ndofs; k++) {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;

            const int i = (dof_iterator)->EquationId();
            (dof_iterator)->GetSolutionStepReactionValue() = -b[i];
        }
    }

    /**
     * @brief Applies the dirichlet conditions. This operation may be very heavy or completely
     * unexpensive depending on the implementation choosen and on how the System Matrix is built.
     * @details For explanation of how it works for a particular implementation the user
     * should refer to the particular Builder And Solver choosen
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPartType& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        const double start_apply = OpenMPUtils::GetCurrentTime();

        std::size_t system_size = A.size1();
        std::vector<ValueType> scaling_factors(system_size, 0.0);

        // TODO ndofs should be size_t here
        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(ndofs)
        for (int k = 0; k<ndofs; k++) {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;
            if(dof_iterator->IsFixed())
                scaling_factors[k] = 0.0;
            else
                scaling_factors[k] = 1.0;
        }

        auto* Avalues = A.value_data().begin();
        auto* Arow_indices = A.index1_data().begin();
        auto* Acol_indices = A.index2_data().begin();

        //detect if there is a line of all zeros and set the diagonal to a 1 if this happens
        #pragma omp parallel for firstprivate(system_size)
        for (int k = 0; k < static_cast<int>(system_size); ++k){
            auto col_begin = Arow_indices[k];
            auto col_end = Arow_indices[k+1];
            bool empty = true;
            for (auto j = col_begin; j < col_end; ++j)
            {
                if(std::abs(Avalues[j]) > 1.0e-13)
                {
                    empty = false;
                    break;
                }
            }

            if(empty == true)
            {
                A(k,k) = 1.0;
                b[k] = 0.0;
            }
        }

        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(system_size); ++k)
        {
            auto col_begin = Arow_indices[k];
            auto col_end = Arow_indices[k+1];
            ValueType k_factor = scaling_factors[k];
            if (k_factor == 0)
            {
                // zero out the whole row, except the diagonal
                for (auto j = col_begin; j < col_end; ++j)
                    if (static_cast<int>(Acol_indices[j]) != k )
                        Avalues[j] = 0.0;

                // zero out the RHS
                b[k] = 0.0;
            }
            else
            {
                // zero out the column which is associated with the zero'ed row
                for (auto j = col_begin; j < col_end; ++j)
                    if(scaling_factors[ Acol_indices[j] ] == 0 )
                        Avalues[j] = 0.0;
            }
        }

        const double stop_apply = OpenMPUtils::GetCurrentTime();

        if (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "Apply Dirichlet time: " << stop_apply - start_apply << std::endl;
        }
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
     */
    void Clear() override
    {
        BaseType::Clear();
        this->mpLinearSystemSolver->Clear();

        mSlaveIds.clear();
        mMasterIds.clear();
        mInactiveSlaveDofs.clear();
        mT.resize(0,0,false);
        mConstantVector.resize(0,false);

        if (this->GetEchoLevel() > 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints Clear Function called" << std::endl;
        }
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model part of the problem to solve
     * @return 0 all ok
     */
    int Check(const ModelPartType& rModelPart) const override
    {
        KRATOS_TRY

        return 0;
        KRATOS_CATCH("");
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedBlockBuilderAndSolverWithConstraints";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    TSystemMatrixType mT;              /// This is matrix containing the global relation for the constraints
    TSystemVectorType mConstantVector; /// This is vector containing the rigid movement of the constraint
    std::vector<IndexType> mSlaveIds;  /// The equation ids of the slaves
    std::vector<IndexType> mMasterIds; /// The equation ids of the master
    std::unordered_set<IndexType> mInactiveSlaveDofs; /// The set containing the inactive slave dofs

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void ConstructMasterSlaveConstraintsStructure(ModelPartType& rModelPart)
    {
        if (rModelPart.MasterSlaveConstraints().size() > 0) {
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // Vector containing the localization in the system of the different terms
            DofsVectorType slave_dof_list, master_dof_list;

            // Constraint initial iterator
            const auto it_const_begin = rModelPart.MasterSlaveConstraints().begin();
            std::vector<std::unordered_set<IndexType>> indices(BaseType::mDofSet.size());

            std::vector<omp_lock_t> lock_array(indices.size());

            for(std::size_t i = 0; i < indices.size(); ++i)
                omp_init_lock(&lock_array[i]);

            #pragma omp parallel firstprivate(slave_dof_list, master_dof_list)
            {
                typename ElementType::EquationIdVectorType slave_ids(3);
                typename ElementType::EquationIdVectorType master_ids(3);
                std::unordered_map<IndexType, std::unordered_set<IndexType>> temp_indices;

                #pragma omp for schedule(guided, 512) nowait
                for (int i_const = 0; i_const < static_cast<int>(rModelPart.MasterSlaveConstraints().size()); ++i_const) {
                    auto it_const = it_const_begin + i_const;

                    // Detect if the constraint is active or not. If the user did not make any choice the constraint
                    // It is active by default
                    bool constraint_is_active = true;
                    if( it_const->IsDefined(ACTIVE) ) {
                        constraint_is_active = it_const->Is(ACTIVE);
                    }

                    if(constraint_is_active) {
                        it_const->EquationIdVector(slave_ids, master_ids, r_current_process_info);

                        // Slave DoFs
                        for (auto &id_i : slave_ids) {
                            temp_indices[id_i].insert(master_ids.begin(), master_ids.end());
                        }
                    }
                }

                // Merging all the temporal indexes
                for (int i = 0; i < static_cast<int>(temp_indices.size()); ++i) {
                    omp_set_lock(&lock_array[i]);
                    indices[i].insert(temp_indices[i].begin(), temp_indices[i].end());
                    omp_unset_lock(&lock_array[i]);
                }
            }

            for(std::size_t i = 0; i < indices.size(); ++i)
                omp_destroy_lock(&lock_array[i]);

            mSlaveIds.clear();
            mMasterIds.clear();
            for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
                if (indices[i].size() == 0) // Master dof!
                    mMasterIds.push_back(i);
                else // Slave dof
                    mSlaveIds.push_back(i);
                indices[i].insert(i); // Ensure that the diagonal is there in T
            }

            // Count the row sizes
            std::size_t nnz = 0;
            for (IndexType i = 0; i < indices.size(); ++i)
                nnz += indices[i].size();

            mT = TSystemMatrixType(indices.size(), indices.size(), nnz);
            mConstantVector.resize(indices.size(), false);

            auto *Tvalues = mT.value_data().begin();
            IndexType *Trow_indices = mT.index1_data().begin();
            IndexType *Tcol_indices = mT.index2_data().begin();

            // Filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
            Trow_indices[0] = 0;
            for (int i = 0; i < static_cast<int>(mT.size1()); i++)
                Trow_indices[i + 1] = Trow_indices[i] + indices[i].size();

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(mT.size1()); ++i) {
                const IndexType row_begin = Trow_indices[i];
                const IndexType row_end = Trow_indices[i + 1];
                IndexType k = row_begin;
                for (auto it = indices[i].begin(); it != indices[i].end(); ++it) {
                    Tcol_indices[k] = *it;
                    Tvalues[k] = 0.0;
                    k++;
                }

                indices[i].clear(); //deallocating the memory

                std::sort(&Tcol_indices[row_begin], &Tcol_indices[row_end]);
            }

            mT.set_filled(indices.size() + 1, nnz);

            Timer::Stop("ConstraintsRelationMatrixStructure");
        }
    }

    void BuildMasterSlaveConstraints(ModelPartType& rModelPart)
    {
        KRATOS_TRY

        TSparseSpace::SetToZero(mT);
        TSparseSpace::SetToZero(mConstantVector);

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Vector containing the localization in the system of the different terms
        DofsVectorType slave_dof_list, master_dof_list;

        // Contributions to the system
        LocalSystemMatrixType transformation_matrix(0, 0);
        LocalSystemVectorType constant_vector(0);

        // Vector containing the localization in the system of the different terms
        typename ElementType::EquationIdVectorType slave_equation_ids, master_equation_ids;

        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

        // We clear the set
        mInactiveSlaveDofs.clear();

        #pragma omp parallel firstprivate(transformation_matrix, constant_vector, slave_equation_ids, master_equation_ids)
        {
            std::unordered_set<IndexType> auxiliar_inactive_slave_dofs;

            #pragma omp for schedule(guided, 512)
            for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
                auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;

                // Detect if the constraint is active or not. If the user did not make any choice the constraint
                // It is active by default
                bool constraint_is_active = true;
                if (it_const->IsDefined(ACTIVE))
                    constraint_is_active = it_const->Is(ACTIVE);

                if (constraint_is_active) {
                    it_const->CalculateLocalSystem(transformation_matrix, constant_vector, r_current_process_info);
                    it_const->EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);

                    for (IndexType i = 0; i < slave_equation_ids.size(); ++i) {
                        const IndexType i_global = slave_equation_ids[i];

                        // Assemble matrix row
                        AssembleRowContribution(mT, transformation_matrix, i_global, i, master_equation_ids);

                        // Assemble constant vector
                        const auto constant_value = constant_vector[i];
                        auto& r_value = mConstantVector[i_global];
                        if constexpr (std::is_arithmetic<TDataType>::value)
                        {
                            #pragma omp atomic
                            r_value += constant_value;
                        }
                        else
                        {
                            #pragma omp critical
                            {
                                r_value += constant_value;
                            }
                        }
                    }
                } else { // Taking into account inactive constraints
                    it_const->EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);
                    auxiliar_inactive_slave_dofs.insert(slave_equation_ids.begin(), slave_equation_ids.end());
                }
            }

            // We merge all the sets in one thread
            #pragma omp critical
            {
                mInactiveSlaveDofs.insert(auxiliar_inactive_slave_dofs.begin(), auxiliar_inactive_slave_dofs.end());
            }
        }

        // Setting the master dofs into the T and C system
        for (auto eq_id : mMasterIds) {
            mConstantVector[eq_id] = 0.0;
            mT(eq_id, eq_id) = 1.0;
        }

        // Setting inactive slave dofs in the T and C system
        for (auto eq_id : mInactiveSlaveDofs) {
            mConstantVector[eq_id] = 0.0;
            mT(eq_id, eq_id) = 1.0;
        }

        // fixing the constraint transformation matrix if there's a dof that's a slave in one constraint and a master in another constraint.
        // TODO to be reviewed (code from Yaman)
        for (int k = 0; k < static_cast<int>(mSlaveIds.size()); ++k)
        {
            const IndexType slave_equation_id = mSlaveIds[k];
            for (int i = 0; i < static_cast<int>(mT.size1()); ++i)
            {
                TDataType tmp = mT(i, slave_equation_id);
                if (std::abs(tmp) > 1.0e-14)
                {
                    std::cout << "ATTENTION! ResidualBasedBlockBuilderAndSolverWithConstraints. constraint slave is used as master for another constraint!" << std::endl;
                    for (int j = 0; j < static_cast<int>(mT.size2()); ++j)
                    {
                        TDataType v1 = mT(i, slave_equation_id);
                        TDataType v2 = mT(slave_equation_id, j);
                        mT(i, j) += v1*v2;
                    }
                    mT(i, slave_equation_id) = 0.0;
                }
            }
        }

        KRATOS_CATCH("")
    }

    void ApplyRHSConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPartType& rModelPart,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            double time_begin = OpenMPUtils::GetCurrentTime();
            double time_end, time_1, time_2, time_3, time_4;

            BuildMasterSlaveConstraints(rModelPart);

            time_end = OpenMPUtils::GetCurrentTime();
            time_1 = time_end - time_begin;
            time_begin = time_end;

            // We compute the transposed matrix of the global relation matrix
            TSystemMatrixType T_transpose_matrix(mT.size2(), mT.size1());
            SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, mT, 1.0);

            time_end = OpenMPUtils::GetCurrentTime();
            time_2 = time_end - time_begin;
            time_begin = time_end;

            TSystemVectorType b_modified(rb.size());
            TSparseSpace::Mult(T_transpose_matrix, rb, b_modified);
            TSparseSpace::Copy(b_modified, rb);
            b_modified.resize(0, false); //free memory

            time_end = OpenMPUtils::GetCurrentTime();
            time_3 = time_end - time_begin;
            time_begin = time_end;

            // Apply diagonal values on slaves
            // #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(mSlaveIds.size()); ++i) {
                const IndexType slave_equation_id = mSlaveIds[i];
                if (mInactiveSlaveDofs.find(slave_equation_id) == mInactiveSlaveDofs.end()) {
                    rb[slave_equation_id] = 0.0;
                }
            }

            time_end = OpenMPUtils::GetCurrentTime();
            time_4 = time_end - time_begin;

            if (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)
            {
                std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "Apply RHS Constraints time: ("
                          << time_1 << ", " << time_2 << ", " << time_3 << ", " << time_4 << ")"
                          << ", total = " << time_1 + time_2 + time_3 + time_4
                          << std::endl;
            }
        }

        KRATOS_CATCH("")
    }

    void ApplyConstraints(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType &rA,
        TSystemVectorType &rDx,
        TSystemVectorType &rb,
        ModelPartType &rModelPart)
    {
        KRATOS_TRY

        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            const double start_apply = OpenMPUtils::GetCurrentTime();

            BuildMasterSlaveConstraints(rModelPart);

            // We compute the transposed matrix of the global relation matrix
            TSystemMatrixType T_transpose_matrix(mT.size2(), mT.size1());
            SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, mT, 1.0);

            TSystemVectorType b_modified(rb.size());
            TSparseSpace::Mult(T_transpose_matrix, rb, b_modified);
            TSparseSpace::Copy(b_modified, rb);
            b_modified.resize(0, false); //free memory

            TSystemMatrixType auxiliar_A_matrix(mT.size2(), rA.size2());
            SparseMatrixMultiplicationUtility::MatrixMultiplication(T_transpose_matrix, rA, auxiliar_A_matrix); //auxiliar = T_transpose * rA
            T_transpose_matrix.resize(0, 0, false);                                                             //free memory

            SparseMatrixMultiplicationUtility::MatrixMultiplication(auxiliar_A_matrix, mT, rA); //A = auxilar * T   NOTE: here we are overwriting the old A matrix!
            auxiliar_A_matrix.resize(0, 0, false);                                              //free memory

            ValueType max_diag = 0.0;
            for(IndexType i = 0; i < rA.size1(); ++i) {
                TDataType tmp = rA(i, i);
                max_diag = std::max(std::abs(tmp), max_diag);
            }

            // Apply diagonal values on slaves
            // #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(mSlaveIds.size()); ++i) {
                const IndexType slave_equation_id = mSlaveIds[i];
                if (mInactiveSlaveDofs.find(slave_equation_id) == mInactiveSlaveDofs.end()) {
                    rA(slave_equation_id, slave_equation_id) = max_diag;
                    rb[slave_equation_id] = 0.0;
                }
            }

            const double stop_apply = OpenMPUtils::GetCurrentTime();

            if (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)
            {
                std::cout << "ResidualBasedBlockBuilderAndSolverWithConstraints: " << "Apply Constraints time: " << stop_apply - start_apply << std::endl;
            }
        }

        KRATOS_CATCH("")
    }

    virtual void ConstructMatrixStructure(
        TSystemMatrixType& A,
        const ModelPartType& rModelPart) const
    {
        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        auto el_begin = rModelPart.ElementsBegin();
        auto cond_begin = rModelPart.ConditionsBegin();

        const SizeType equation_size = BaseType::mEquationSystemSize;

        std::vector< omp_lock_t > lock_array(equation_size);

        for(IndexType i = 0; i < equation_size; ++i)
            omp_init_lock(&lock_array[i]);

        std::vector<std::unordered_set<IndexType> > indices(equation_size);
        #pragma omp parallel for firstprivate(equation_size)
        for (int iii = 0; iii < static_cast<int>(equation_size); iii++) {
            indices[iii].reserve(40);
        }

        typename ElementType::EquationIdVectorType ids(3, 0);

        #pragma omp parallel for firstprivate(nelements, ids)
        for (int iii=0; iii<nelements; iii++) {
            auto i_element = el_begin + iii;
            i_element->EquationIdVector(ids, CurrentProcessInfo);
            for (std::size_t i = 0; i < ids.size(); i++) {
                omp_set_lock(&lock_array[ids[i]]);
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
                omp_unset_lock(&lock_array[ids[i]]);
            }
        }

        #pragma omp parallel for firstprivate(nconditions, ids)
        for (int iii = 0; iii<nconditions; iii++) {
            auto i_condition = cond_begin + iii;
            i_condition->EquationIdVector(ids, CurrentProcessInfo);
            for (std::size_t i = 0; i < ids.size(); i++) {
                omp_set_lock(&lock_array[ids[i]]);
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
                omp_unset_lock(&lock_array[ids[i]]);
            }
        }

        //destroy locks
        for(std::size_t i = 0; i < equation_size; ++i)
            omp_destroy_lock(&lock_array[i]);

        //count the row sizes
        unsigned int nnz = 0;
        for (unsigned int i = 0; i < indices.size(); i++) {
            nnz += indices[i].size();
        }

        A = CompressedMatrixType(indices.size(), indices.size(), nnz);

        auto* Avalues = A.value_data().begin();
        auto* Arow_indices = A.index1_data().begin();
        auto* Acol_indices = A.index2_data().begin();

        //filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
        Arow_indices[0] = 0;
        for (int i = 0; i < static_cast<int>(A.size1()); i++) {
            Arow_indices[i+1] = Arow_indices[i] + indices[i].size();
        }

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(A.size1()); i++) {
            const unsigned int row_begin = Arow_indices[i];
            const unsigned int row_end = Arow_indices[i+1];
            unsigned int k = row_begin;
            for (auto it = indices[i].begin(); it != indices[i].end(); it++) {
                Acol_indices[k] = *it;
                Avalues[k] = 0.0;
                k++;
            }

            indices[i].clear(); //deallocating the memory

            std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);

        }

        A.set_filled(indices.size()+1, nnz);

        Timer::Stop("MatrixStructure");
    }

    void Assemble(
        TSystemMatrixType& A,
        TSystemVectorType& b,
        const LocalSystemMatrixType& LHS_Contribution,
        const LocalSystemVectorType& RHS_Contribution,
        typename ElementType::EquationIdVectorType& EquationId
    ) const
    {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++) {
            unsigned int i_global = EquationId[i_local];

            auto& r_a = b[i_global];
            const auto& v_a = RHS_Contribution(i_local);
            if constexpr (std::is_arithmetic<TDataType>::value)
            {
                #pragma omp atomic
                r_a += v_a;
            }
            else
            {
                #pragma omp critical
                {
                    r_a += v_a;
                }
            }

            AssembleRowContribution(A, LHS_Contribution, i_global, i_local, EquationId);
        }
    }


    //**************************************************************************

    void AssembleRHS(
        TSystemVectorType& b,
        LocalSystemVectorType& RHS_Contribution,
        typename ElementType::EquationIdVectorType& EquationId
    ) const
    {
        unsigned int local_size = RHS_Contribution.size();

        for (unsigned int i_local = 0; i_local < local_size; i_local++) {
            unsigned int i_global = EquationId[i_local];

            // ASSEMBLING THE SYSTEM VECTOR
            auto& b_value = b[i_global];
            const auto& rhs_value = RHS_Contribution[i_local];

            if constexpr (std::is_arithmetic<TDataType>::value)
            {
                #pragma omp atomic
                b_value += rhs_value;
            }
            else
            {
                #pragma omp critical
                {
                    b_value += rhs_value;
                }
            }
        }
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate) const
    {
        std::vector<std::size_t>::iterator i = v.begin();
        std::vector<std::size_t>::iterator endit = v.end();
        while (i != endit && (*i) != candidate) {
            i++;
        }
        if (i == endit) {
            v.push_back(candidate);
        }

    }

    void BuildRHSNoDirichlet(
        typename TSchemeType::Pointer pScheme,
        ModelPartType& rModelPart,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        //Getting the Elements
        ElementsContainerType& pElements = rModelPart.Elements();

        //getting the array of the conditions
        ConditionsContainerType& ConditionsArray = rModelPart.Conditions();

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        typename ElementType::EquationIdVectorType EquationId;

        // assemble all elements
        //for (typename ElementsContainerType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)

        const int nelements = static_cast<int>(pElements.size());
        #pragma omp parallel firstprivate(nelements, RHS_Contribution, EquationId)
        {
            #pragma omp for schedule(guided, 512) nowait
            for (int i=0; i<nelements; i++) {
                auto it = pElements.begin() + i;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if( (it)->IsDefined(ACTIVE) )
                    element_is_active = (it)->Is(ACTIVE);

                if(element_is_active) {
                    //calculate elemental Right Hand Side Contribution
                    pScheme->CalculateRHSContribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }

            LHS_Contribution.resize(0, 0, false);
            RHS_Contribution.resize(0, false);

            // assemble all conditions
            const int nconditions = static_cast<int>(ConditionsArray.size());
            #pragma omp for schedule(guided, 512)
            for (int i = 0; i<nconditions; i++) {
                auto it = ConditionsArray.begin() + i;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool condition_is_active = true;
                if( (it)->IsDefined(ACTIVE) )
                    condition_is_active = (it)->Is(ACTIVE);

                if(condition_is_active) {
                    //calculate elemental contribution
                    pScheme->CalculateRHSContribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }
        }

        KRATOS_CATCH("")
    }

    //******************************************************************************************
    //******************************************************************************************

    inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions) const
    {
        partitions.resize(number_of_threads + 1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for (unsigned int i = 1; i < number_of_threads; i++) {
            partitions[i] = partitions[i - 1] + partition_size;
        }
    }

    inline void AssembleRowContribution(TSystemMatrixType& A, const LocalSystemMatrixType& Alocal,
            const unsigned int i, const unsigned int i_local, typename ElementType::EquationIdVectorType& EquationId) const
    {
        if (EquationId.size() == 0)
            return;

        auto* values_vector = A.value_data().begin();
        auto* index1_vector = A.index1_data().begin();
        auto* index2_vector = A.index2_data().begin();

        IndexType left_limit = index1_vector[i];
//    size_t right_limit = index1_vector[i+1];

        //find the first entry
        IndexType last_pos = ForwardFind(EquationId[0],left_limit,index2_vector);
        IndexType last_found = EquationId[0];

        auto& r_a = values_vector[last_pos];
        const auto& v_a = Alocal(i_local,0);
        if constexpr (std::is_arithmetic<TDataType>::value)
        {
            #pragma omp atomic
            r_a += v_a;
        }
        else
        {
            #pragma omp critical
            {
                r_a += v_a;
            }
        }

        //now find all of the other entries
        IndexType pos = 0;
        for (unsigned int j=1; j<EquationId.size(); j++) {
            unsigned int id_to_find = EquationId[j];
            if(id_to_find > last_found) {
                pos = ForwardFind(id_to_find,last_pos+1,index2_vector);
            } else if(id_to_find < last_found) {
                pos = BackwardFind(id_to_find,last_pos-1,index2_vector);
            } else {
                pos = last_pos;
            }

            auto& r = values_vector[pos];
            const auto& v = Alocal(i_local,j);
            if constexpr (std::is_arithmetic<TDataType>::value)
            {
                #pragma omp atomic
                r += v;
            }
            else
            {
                #pragma omp critical
                {
                    r += v;
                }
            }

            last_found = id_to_find;
            last_pos = pos;
        }
    }

    inline unsigned int ForwardFind(const unsigned int id_to_find,
                                    const unsigned int start,
                                    const size_t* index_vector) const
    {
        unsigned int pos = start;
        while(id_to_find != index_vector[pos]) pos++;
        return pos;
    }

    inline unsigned int BackwardFind(const unsigned int id_to_find,
                                     const unsigned int start,
                                     const size_t* index_vector) const
    {
        unsigned int pos = start;
        while(id_to_find != index_vector[pos]) pos--;
        return pos;
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class ResidualBasedBlockBuilderAndSolverWithConstraints */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

#undef ENABLE_PARALLEL_SETUP_DOFSET

#endif /* KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS  defined */
