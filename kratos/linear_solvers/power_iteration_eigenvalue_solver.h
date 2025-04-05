//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-15 11:11:35 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_POWER_ITERATION_EIGENVALUE_SOLVERR_H_INCLUDED )
#define  KRATOS_POWER_ITERATION_EIGENVALUE_SOLVERR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <numeric>
#include <vector>


// External includes


// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"


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

/// Short class definition.
/** Detail class definition.
*/
template<class TSparseSpaceType, class TDenseSpaceType, class TLinearSolverType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class PowerIterationEigenvalueSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PowerIterationEigenvalueSolver
    KRATOS_CLASS_POINTER_DEFINITION(PowerIterationEigenvalueSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::DenseMatrixType DenseMatrixType;

    typedef typename BaseType::DenseVectorType DenseVectorType;

    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::IndexType IndexType;

    typedef typename BaseType::DataType DataType;

    typedef typename BaseType::ValueType ValueType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PowerIterationEigenvalueSolver() {}

    PowerIterationEigenvalueSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber,
                                   unsigned int NewRequiredEigenvalueNumber, typename TLinearSolverType::Pointer pLinearSolver)
        : BaseType(NewMaxTolerance, NewMaxIterationsNumber), mRequiredEigenvalueNumber(NewRequiredEigenvalueNumber), mpLinearSolver(pLinearSolver) {}

    /*       PowerIterationEigenvalueSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) :  */
    /*       BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){} */

    /// Copy constructor.
    PowerIterationEigenvalueSolver(const PowerIterationEigenvalueSolver& Other) : BaseType(Other) {}


    /// Destructor.
    virtual ~PowerIterationEigenvalueSolver() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    PowerIterationEigenvalueSolver& operator=(const PowerIterationEigenvalueSolver& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    static void RandomInitialize(DenseVectorType& R)
    {
        for(SizeType i = 0 ; i < R.size() ; i++)
            R[i] = 1.00; //rand();

        R /= norm_2(R);
    }


    // The power iteration algorithm
    void Solve(SparseMatrixType& K,
               SparseMatrixType& M,
               DenseVectorType& Eigenvalues,
               DenseMatrixType& Eigenvectors)
    {

        using boost::numeric::ublas::trans;

        SizeType size = K.size1();
        SizeType max_iteration = BaseType::GetMaxIterationsNumber();
        ValueType tolerance = BaseType::GetTolerance();

        VectorType x = ZeroVector(size);
        VectorType y = ZeroVector(size);

        RandomInitialize(y);

        if(Eigenvalues.size() < 1)
            Eigenvalues.resize(1,0.00);


        // Starting with first step
        DataType beta = 0.00;
        DataType ro = 0.00;
        DataType old_ro = Eigenvalues[0];
        std::cout << "iteration    beta \t\t ro \t\t convergence norm" << std::endl;
        for(IndexType i = 0 ; i < max_iteration ; i++)
        {
            //K*x = y
            mpLinearSolver->Solve(K,x,y);

            ro = inner_prod(y,x);

            //y = M*x
            noalias(y) = prod(M,x);

            beta = inner_prod(x, y);
            if constexpr (std::is_same<DataType, ValueType>::value)
            {
                if(beta <= 0.00)
                    KRATOS_THROW_ERROR(std::invalid_argument, "M is not Positive-definite", "");
            }

            ro = ro / beta;
            beta = std::sqrt(beta);

            DataType inverse_of_beta = 1.00 / beta;

            y *= inverse_of_beta;

            if(std::abs(ro) == 0.00)
                KRATOS_THROW_ERROR(std::runtime_error, "Perpendicular eigenvector to M", "");

            ValueType convergence_norm = std::abs((ro - old_ro) / ro);

            std::cout << i << " \t " << beta << " \t " << ro << " \t " << convergence_norm << std::endl;
            //std::cout << "i = " << i << ": beta = " << beta << ", ro = " << ro << ", convergence norm = " << convergence_norm << std::endl;

            if(convergence_norm < tolerance)
                break;

            old_ro = ro;



        }

        KRATOS_WATCH(ro);
//KRATOS_WATCH(y);

        Eigenvalues[0] = ro;

        if((Eigenvectors.size1() < 1) || (Eigenvectors.size2() < size))
            Eigenvectors.resize(1,size);

        for(SizeType i = 0 ; i < size ; i++)
            Eigenvectors(0,i) = y[i];
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
        std::stringstream buffer;
        buffer << "Power iteration eigenvalue solver with " << BaseType::GetPreconditioner()->Info();
        return  buffer.str();
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


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


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


    unsigned int mRequiredEigenvalueNumber;

    typename TLinearSolverType::Pointer mpLinearSolver;

    std::vector<DenseVectorType> mQVector;
    std::vector<DenseVectorType> mPVector;
    std::vector<DenseVectorType> mRVector;

    ///@}
    ///@name Private Operators
    ///@{


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

}; // Class PowerIterationEigenvalueSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_POWER_ITERATION_EIGENVALUE_SOLVERR_H_INCLUDED defined
