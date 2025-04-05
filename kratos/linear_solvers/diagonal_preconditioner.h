//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_DIAGONAL_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_DIAGONAL_PRECONDITIONER_H_INCLUDED



// System includes


// External includes

// Project includes
#include "includes/define.h"


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

///@name  Preconditioners
///@{

/// DiagonalPreconditioner class.
/** DiagonalPreconditioner for linesr system solvers.
 */
template<class TSparseSpaceType, class TDenseSpaceType>
class DiagonalPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of DiagonalPreconditioner
    KRATOS_CLASS_POINTER_DEFINITION(DiagonalPreconditioner);

    typedef  Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef typename TSparseSpaceType::DataType DataType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DiagonalPreconditioner() {}

    /// Copy constructor.
    DiagonalPreconditioner(const DiagonalPreconditioner& Other)
        : BaseType(Other), mDiagonal(Other.mDiagonal), mTemp(Other.mTemp) {}

    /// Destructor.
    virtual ~DiagonalPreconditioner() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    DiagonalPreconditioner& operator=(const DiagonalPreconditioner& Other)
    {
        BaseType::operator=(Other);
        mDiagonal = Other.mDiagonal;
        mTemp = Other.mTemp;
        return *this;
    }


    ///@}
    ///@name Operations
    ///@{


    /** DiagonalPreconditioner Initialize
    Initialize preconditioner for linear system rA*rX=rB
    @param rA  system matrix.
    @param rX Unknows vector
    @param rB Right side linear system of equations.
    */
    void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        mDiagonal.resize(int(TSparseSpaceType::Size(rX)));
        mTemp.resize(int(TSparseSpaceType::Size(rX)));

        int i;

        const DataType zero = DataType();
// 	  const DataType one = DataType(1.00);

        #pragma omp parallel for private(i)
        for(i = 0 ; i < int(rA.size1()) ; ++i)
        {
            DataType diag_Aii = rA(i,i);
            if(diag_Aii != zero)
                mDiagonal[i] = 1.00 / sqrt(fabs(diag_Aii));
            else
                KRATOS_THROW_ERROR(std::logic_error,"zero found in the diagonal. Diagonal preconditioner can not be used","");
        }
// 	      mDiagonal[i] = one;

        /* 	  std::cout << "mDiagonal : " << mDiagonal << std::endl; */

        /* for(i = 0 ; i < rA.RowsNumber() ; ++i)
           for(j = 0 ; j < rA.RowsNumber() ; ++j) */


    }

    void Initialize(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
        BaseType::Initialize(rA, rX, rB);
    }

    void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
    {
        int i;
        #pragma omp parallel for private(i)
        for(i = 0 ; i < int(TSparseSpaceType::Size(rX)) ; ++i)
            mTemp[i] = rX[i] * mDiagonal[i];
        TSparseSpaceType::Mult(rA,mTemp, rY);
        ApplyLeft(rY);
    }

    void TransposeMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
    {
        int i;
        #pragma omp parallel for private(i)
        for(i = 0 ; i < int(TSparseSpaceType::Size(rX)) ; ++i)
            mTemp[i] = rX[i] * mDiagonal[i];
        TSparseSpaceType::TransposeMult(rA,mTemp, rY);
        ApplyRight(rY);
    }

    VectorType& ApplyLeft(VectorType& rX)
    {
        int i;
        #pragma omp parallel for private(i)
        for(i = 0 ; i < int(TSparseSpaceType::Size(rX)) ; ++i)
            rX[i] *= mDiagonal[i];

        return rX;
    }

    VectorType& ApplyRight(VectorType& rX)
    {
        int i;
        #pragma omp parallel for private(i)
        for(i = 0 ; i < int(TSparseSpaceType::Size(rX)) ; ++i)
            rX[i] *= mDiagonal[i];

        return rX;
    }

    /** DiagonalPreconditioner transpose solver.
    Solving tranpose preconditioner system M^T*x=y, where m^T means transpose.
    @param rMatrix   DiagonalPreconditioner system matrix.
    @param rXVector  Unknows of preconditioner suystem
    @param rYVector  Right side of preconditioner system.
    */
    VectorType& ApplyTransposeLeft(VectorType& rX)
    {
        int i;
        #pragma omp parallel for private(i)
        for(i = 0 ; i < int(TSparseSpaceType::Size(rX)) ; ++i)
            rX[i] *= mDiagonal[i];

        return rX;
    }

    VectorType& ApplyTransposeRight(VectorType& rX)
    {
        int i;
        #pragma omp parallel for private(i)
        for(i = 0 ; i < int(TSparseSpaceType::Size(rX)) ; ++i)
            rX[i] *= mDiagonal[i];

        return rX;
    }

    VectorType& ApplyInverseRight(VectorType& rX)
    {
// 	  const DataType zero = DataType();

        int i;
        #pragma omp parallel for private(i)
        for(i = 0 ; i < int(TSparseSpaceType::Size(rX)) ; ++i)
            /*	    if(mDiagonal[i] != zero)*/
            rX[i] /= mDiagonal[i];

        return rX;
    }

    VectorType& Finalize(VectorType& rX)
    {
        int i;
        #pragma omp parallel for private(i)
        for(i = 0 ; i < int(TSparseSpaceType::Size(rX)) ; ++i)
            rX[i] *= mDiagonal[i];

        return rX;
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

    /// Return information about this object.
    std::string  Info() const override
    {
        return "DiagonalPreconditioner";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& OStream) const override
    {
        OStream << "Diagonal preconditioner";
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

    VectorType mDiagonal;

    VectorType mTemp;


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

}; // Class DiagonalPreconditioner

///@}

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_DIAGONAL_PRECONDITIONER_H_INCLUDED  defined
