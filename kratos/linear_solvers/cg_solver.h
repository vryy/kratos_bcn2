//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_CG_SOLVER_H_INCLUDED )
#define  KRATOS_CG_SOLVER_H_INCLUDED



// System includes
#include <string>
#include <iostream>


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
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class CGSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CGSolver
    KRATOS_CLASS_POINTER_DEFINITION(CGSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename BaseType::DataType DataType;

    typedef typename BaseType::ValueType ValueType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::DenseMatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CGSolver() {}

    CGSolver(ValueType NewMaxTolerance) : BaseType(NewMaxTolerance) {}

    CGSolver(ValueType NewMaxTolerance, unsigned int NewMaxIterationsNumber) : BaseType(NewMaxTolerance, NewMaxIterationsNumber) {}

    CGSolver(ValueType NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) :
        BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner) {}

    /// Copy constructor.
    CGSolver(const CGSolver& Other) : BaseType(Other) {}


    /// Destructor.
    virtual ~CGSolver() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    CGSolver& operator=(const CGSolver& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /** Normal solve method.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

// 	  GetTimeTable()->Start(Info());

        BaseType::GetPreconditioner()->Initialize(rA,rX,rB);
        BaseType::GetPreconditioner()->ApplyInverseRight(rX);
        BaseType::GetPreconditioner()->ApplyLeft(rB);

        bool is_solved = IterativeSolve(rA,rX,rB);

        BaseType::GetPreconditioner()->Finalize(rX);

// 	  GetTimeTable()->Stop(Info());

        return is_solved;
    }


    /** Multi solve method for solving a set of linear systems with same coefficient matrix.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
// 	  GetTimeTable()->Start(Info());

        BaseType::GetPreconditioner()->Initialize(rA,rX,rB);

        bool is_solved = true;
        VectorType x(TDenseSpaceType::Size1(rX));
        VectorType b(TDenseSpaceType::Size1(rB));
        for(unsigned int i = 0 ; i < TDenseSpaceType::Size2(rX) ; i++)
        {
            TDenseSpaceType::GetColumn(i,rX, x);
            TDenseSpaceType::GetColumn(i,rB, b);

            BaseType::GetPreconditioner()->ApplyInverseRight(x);
            BaseType::GetPreconditioner()->ApplyLeft(b);

            is_solved &= IterativeSolve(rA,x,b);

            BaseType::GetPreconditioner()->Finalize(x);
        }

// 	  GetTimeTable()->Stop(Info());

        return is_solved;
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
        buffer << "Conjugate gradient linear solver with " << BaseType::GetPreconditioner()->Info();
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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    bool IterativeSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        const int size = TSparseSpaceType::Size(rX);

        BaseType::mIterationsNumber = 0;

        VectorType r(size);

        this->PreconditionedMult(rA,rX,r);
        TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);

        BaseType::mBNorm = std::abs(TSparseSpaceType::TwoNorm(rB));

        VectorType p(r);
        VectorType q(size);

        DataType roh0 = TSparseSpaceType::Dot(r, r);
        DataType roh1 = roh0;
        DataType beta = 0;

        if(std::abs(roh0) < 1.0e-30) //modification by Riccardo
//	if(roh0 == 0.00)
            return false;

        do
        {
            this->PreconditionedMult(rA,p,q);

            DataType pq = TSparseSpaceType::Dot(p,q);

            //if(pq == 0.00)
            if(std::abs(pq) <= 1.0e-30)
                break;

            DataType alpha = roh0 / pq;

            TSparseSpaceType::ScaleAndAdd(alpha, p, 1.00, rX);
            TSparseSpaceType::ScaleAndAdd(-alpha, q, 1.00, r);

            roh1 = TSparseSpaceType::Dot(r,r);

            beta = (roh1 / roh0);
            TSparseSpaceType::ScaleAndAdd(1.00, r, beta, p);

            roh0 = roh1;

            BaseType::mResidualNorm = std::abs(std::sqrt(roh1));
            BaseType::mIterationsNumber++;

            if (this->GetEchoLevel() > 0)
            {
                std::cout << "CGSolver iteration #" << BaseType::mIterationsNumber
                          << ", normr = " << sqrt(roh1) << ", tol = " << this->GetTolerance()
                          << std::endl;
            }
        }
        while(BaseType::IterationNeeded() && (std::abs(roh0) > 1.0e-30)/*(roh0 != 0.00)*/);

        return BaseType::IsConverged();
    }

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

}; // Class CGSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_CG_SOLVER_H_INCLUDED  defined
