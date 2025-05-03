//*********************************************************************************
//*   Project Name:        Kratos
//*   Last Modified by:    $Author: rrossi $
//*   Date:                $Date: 2007-03-06 10:30:33 $
//*   Revision:            $Revision: 1.2 $
//*********************************************************************************


#if !defined(KRATOS_TFQMR_SOLVER_H_INCLUDED )
#define  KRATOS_TFQMR_SOLVER_H_INCLUDED



// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include <math.h>


namespace Kratos
{

template<class TSparseSpaceType,
         class TDenseSpaceType,
         class TModelPartType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType, TModelPartType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class TFQMRSolver : public IterativeSolver<TSparseSpaceType,
    TDenseSpaceType,
    TModelPartType,
    TPreconditionerType,
    TReordererType>
{
public:

    /// Counted pointer of TFQMRSolver
    KRATOS_CLASS_POINTER_DEFINITION(TFQMRSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TPreconditionerType, TReordererType > BaseType;

    typedef typename BaseType::DataType DataType;

    typedef typename BaseType::ValueType ValueType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::DenseMatrixType DenseMatrixType;

    static constexpr auto zero = DataType();
    static constexpr auto one = static_cast<DataType>(1.0);

    /// Default constructor.
    TFQMRSolver() {}


    TFQMRSolver(ValueType NewTolerance) : BaseType(NewTolerance) {}


    TFQMRSolver(ValueType NewTolerance,
                unsigned int NewMaxIterationsNumber) : BaseType(NewTolerance, NewMaxIterationsNumber) {}


    TFQMRSolver(ValueType NewMaxTolerance,
                unsigned int NewMaxIterationsNumber,
                typename TPreconditionerType::Pointer pNewPreconditioner) :
        BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner) {}

    /// Copy constructor.
    TFQMRSolver(const TFQMRSolver& Other):BaseType(Other) {}


    /// Destructor.
    ~TFQMRSolver() override {}


    /** Normal solve method.
        Solves the linear system Ax=b and puts the result on SystemVector& rX.
        rX is also th initial guess for iterative methods.
        @param rA. System matrix
        @param rX. Solution vector. it's also the initial
        guess for iterative linear solvers.
        @param rB. Right hand side vector*/
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;


        /*            GetTimeTable()->Start(Info()); */


        BaseType::GetPreconditioner()->Initialize(rA,rX,rB);
        BaseType::GetPreconditioner()->ApplyInverseRight(rX);
        BaseType::GetPreconditioner()->ApplyLeft(rB);


        bool is_solved = IterativeSolve(rA,rX,rB);


        BaseType::GetPreconditioner()->Finalize(rX);


        /*            GetTimeTable()->Stop(Info()); */


        return is_solved;
    }

    /** Multi solve method for solving a set of linear systems with same coefficient matrix.
        Solves the linear system Ax=b and puts the result on SystemVector& rX.
        rX is also th initial guess for iterative methods.
        @param rA. System matrix
        @param rX. Solution vector. it's also the initial
        guess for iterative linear solvers.
        @param rB. Right hand side vector*/
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        /*            GetTimeTable()->Start(Info()); */


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


        /*            GetTimeTable()->Stop(Info()); */


        return is_solved;
    }

    /// Return information about this object.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Tranpose-free QMR linear solver with " << BaseType::GetPreconditioner()->Info();
        return  buffer.str();
    }

    /// Print information about this object.
    void  PrintInfo(std::ostream& OStream) const override
    {
        OStream << "Tranpose-free QMR linear solver with ";
        BaseType::GetPreconditioner()->PrintInfo(OStream);
    }


    /// Print object's data.
    void  PrintData(std::ostream& OStream) const override
    {
        BaseType::PrintData(OStream);
    }

private:

    bool IterativeSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        const int size = TSparseSpaceType::Size(rX);


        int info_count    = 0;
        BaseType::mIterationsNumber = 0;


        BaseType::mBNorm        = std::abs(TSparseSpaceType::TwoNorm(rB));


        ValueType errtol = BaseType::mBNorm*BaseType::GetTolerance();


        VectorType r(rB);
        VectorType w(r);


        VectorType y1(r);
        VectorType y2(size);
        for(typename VectorType::iterator ity2 = y2.begin(); ity2 != y2.end(); ity2++) *ity2 = zero;


        VectorType d(size);
        for(typename VectorType::iterator itd = d.begin(); itd != d.end(); itd++) *itd = zero;


        VectorType v(size);
        this->PreconditionedMult(rA, y1, v);


        VectorType u1(v);
        VectorType u2(size);
        for(typename VectorType::iterator itu2 = u2.begin(); itu2 != u2.end(); itu2++) *itu2 = zero;


        DataType theta = zero;
        DataType eta   = zero;
        DataType sigma = zero;
        DataType alpha = zero;
        DataType c     = zero;


        DataType tau   = TSparseSpaceType::TwoNorm(r);
        DataType rho   = tau*tau;
        DataType rhon  = zero;
        DataType beta  = zero;


        int m = 0;



        while(BaseType::mIterationsNumber < BaseType::GetMaxIterationsNumber() )
        {


            BaseType::mIterationsNumber++;


            sigma = TSparseSpaceType::Dot(r,v);

            if (std::abs(sigma) == zero) break;

            alpha = rho/sigma;



            ////////////////even iteration///////////////////


            m = 2*BaseType::mIterationsNumber-1;


            //w=w-alpha*u1;
            TSparseSpaceType::ScaleAndAdd(-alpha, u1, one, w);


            //d=y1+(theta*theta*eta/alpha)*d;
            TSparseSpaceType::ScaleAndAdd(one, y1, (theta*theta*eta/alpha), d);

            theta = TSparseSpaceType::TwoNorm(w)/tau;
            c     = 1 / std::sqrt(1+theta*theta);
            tau   = tau*theta*c;
            eta   = c*c*alpha;


            //x=x+eta*d;
            TSparseSpaceType::ScaleAndAdd(eta, d, one, rX);

            BaseType::mResidualNorm = std::abs(tau)*std::sqrt((ValueType)(m+1));

            if (BaseType::mResidualNorm <= errtol) break;



            ////////////////odd iteration///////////////////


            m = 2*BaseType::mIterationsNumber;


            //y2=y1-alpha*v;
            TSparseSpaceType::ScaleAndAdd(one, y1, -alpha, v, y2);


            //u2=A*y2;
            this->PreconditionedMult(rA,y2,u2);

            //w=w-alpha*u2;
            TSparseSpaceType::ScaleAndAdd(-alpha, u2, one, w);

            //d=y2+(theta*theta*eta/alpha)*d;
            TSparseSpaceType::ScaleAndAdd(one, y2, (theta*theta*eta/alpha), d);


            theta = TSparseSpaceType::TwoNorm(w)/tau;
            c     = one/std::sqrt(1+theta*theta);
            tau   = tau*theta*c;
            eta   = c*c*alpha;


            //x=x+eta*d;
            TSparseSpaceType::ScaleAndAdd(eta, d, one, rX);


            BaseType::mResidualNorm = std::abs(tau)*std::sqrt((ValueType)(m+1));

            if (BaseType::mResidualNorm <= errtol) break;



            ///////////////////////////////////////////////



            if (rho == zero) break;

            rhon = TSparseSpaceType::Dot(r,w);
            beta = rhon/rho;
            rho  = rhon;


            //y1= w + beta*y2;
            TSparseSpaceType::ScaleAndAdd(one, w, beta, y2, y1);


            //u1 = A*y1;
            this->PreconditionedMult(rA,y1,u1);

            //v=u1+beta*(u2+beta*v);
            TSparseSpaceType::ScaleAndAdd(one, u2, beta, v);
            TSparseSpaceType::ScaleAndAdd(one, u1, beta, v);


            //Print iteration info
            info_count++;
            if (info_count==100)
            {
                std::cout<<"it = "<<BaseType::mIterationsNumber<<"   res = "<<(BaseType::mResidualNorm/BaseType::mBNorm)<<std::endl;
                info_count=0;
            }


        }


        std::cout<<std::endl;
        return BaseType::IsConverged();
    }



    /// Assignment operator.
    TFQMRSolver& operator=(const TFQMRSolver& Other);

}; // ********** Class TFQMRSolver *****************

}  // namespace Kratos.

#endif // KRATOS_TFQMR_SOLVER_H_INCLUDED  defined
