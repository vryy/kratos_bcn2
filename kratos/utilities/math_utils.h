/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


/* *********************************************************
*
*   Last Modified by:    $Author: anonymous $
*   Date:                $Date: 2007-12-13 15:29:11 $
*   Revision:            $Revision: 1.8 $
*
* ***********************************************************/


#if !defined(KRATOS_MATH_UTILS )
#define  KRATOS_MATH_UTILS


/* System includes */


/* External includes */
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>


/* Project includes */
#include "includes/define.h"
#include "containers/array_1d.h"


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */


/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

///Various mathematical utilitiy functions.
/**
* Defines several utility functions
 */
template<class TDataType>
class MathUtils
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef boost::numeric::ublas::matrix<TDataType> MatrixType;

    typedef boost::numeric::ublas::vector<TDataType> VectorType;

    typedef unsigned int IndexType;

    /*@} */
    /**@name Life Cycle
             */
    /*@{ */

    /** Constructor.
             */


    /** Destructor.
             */

    /*@} */
    /**@name Operators
             */
    /*@{ */


    /*@} */
    /**@name Operations */
    /*@{ */

    //***********************************************************************
    //***********************************************************************
    static TDataType Distance(const TDataType& rFirstData, const TDataType& rSecondData)
    {
        return rFirstData.Distance(rSecondData);
    }

    //***********************************************************************
    //***********************************************************************
    template<bool check>// = false>
    static inline TDataType Heron(const TDataType& a, const TDataType& b, const TDataType& c)
    {
        TDataType s = 0.5 * (a + b + c);
        TDataType A2 = s * (s - a) * (s - b) * (s - c);
        if(check == true)
        {
            if(A2 < 0.0)
                KRATOS_THROW_ERROR(std::runtime_error, "The square of area is negative, probably the triangle is in bad shape:", A2)
            else
                return sqrt(A2);
        }
        else
            return sqrt(fabs(A2));
    }

    //***********************************************************************
    //***********************************************************************
    static TDataType Abs(const TDataType& rData)
    {
        return rData > TDataType(0) ? rData : -rData;
    }

    static TDataType Min(const TDataType& rValue1, const TDataType& rValue2)
    {
        return rValue1 > rValue2 ? rValue2 : rValue1;
    }

    static TDataType Max(const TDataType& rValue1, const TDataType& rValue2)
    {
        return rValue1 > rValue2 ? rValue1 : rValue2;
    }

    //***********************************************************************
    /**
        inverts matrices of order 2 and 3
             */
    //***********************************************************************
    static inline void InvertMatrix(const MatrixType& InputMatrix,
                                    MatrixType& InvertedMatrix,
                                    TDataType& InputMatrixDet)
    {
        unsigned int size = InputMatrix.size2();

        if (size==2) InvertMatrix2(InputMatrix,InvertedMatrix,InputMatrixDet);
        else InvertMatrix3(InputMatrix,InvertedMatrix,InputMatrixDet);
    }

    //***********************************************************************
    /**
    Inversion of a 3*3 matrix (no check is performed on matrix size):
        - InputMatrix is the input matrix (unchanged at output)
        - InvertedMatrix is the inverse of the input matrix
        - InputMatrixDet is the determinant of the input matrix
             */
    //***********************************************************************
    static void InvertMatrix3(const MatrixType& InputMatrix,
                              MatrixType& InvertedMatrix,
                              TDataType& InputMatrixDet) //VERIFIED!!!
    {
        KRATOS_TRY
        if(InvertedMatrix.size1() != 3 || InvertedMatrix.size2() != 3)
            InvertedMatrix.resize(3,3);

        //filling the inverted matrix with the algebraic complements
        //first column
        InvertedMatrix(0,0) = InputMatrix(1,1)*InputMatrix(2,2) - InputMatrix(1,2)*InputMatrix(2,1);
        InvertedMatrix(1,0) = -InputMatrix(1,0)*InputMatrix(2,2) + InputMatrix(1,2)*InputMatrix(2,0);
        InvertedMatrix(2,0) = InputMatrix(1,0)*InputMatrix(2,1) - InputMatrix(1,1)*InputMatrix(2,0);

        //second column
        InvertedMatrix(0,1) = -InputMatrix(0,1)*InputMatrix(2,2) + InputMatrix(0,2)*InputMatrix(2,1);
        InvertedMatrix(1,1) = InputMatrix(0,0)*InputMatrix(2,2) - InputMatrix(0,2)*InputMatrix(2,0);
        InvertedMatrix(2,1) = -InputMatrix(0,0)*InputMatrix(2,1) + InputMatrix(0,1)*InputMatrix(2,0);

        //third column
        InvertedMatrix(0,2) = InputMatrix(0,1)*InputMatrix(1,2) - InputMatrix(0,2)*InputMatrix(1,1);
        InvertedMatrix(1,2) = -InputMatrix(0,0)*InputMatrix(1,2) + InputMatrix(0,2)*InputMatrix(1,0);
        InvertedMatrix(2,2) = InputMatrix(0,0)*InputMatrix(1,1) - InputMatrix(0,1)*InputMatrix(1,0);

        //calculation of determinant (of the input matrix)
        InputMatrixDet = InputMatrix(0,0)*InvertedMatrix(0,0) + InputMatrix(0,1)*InvertedMatrix(1,0) + InputMatrix(0,2)*InvertedMatrix(2,0);

        //finalizing the calculation of the inverted matrix
        InvertedMatrix /= InputMatrixDet;
        KRATOS_CATCH("")
    }

    //***********************************************************************
    /**
    Inversion of a 2*2 matrix (no check is performed on matrix size):
        - InputMatrix is the input matrix (unchanged at output)
        - InvertedMatrix is the inverse of the input matrix
        - InputMatrixDet is the determinant of the input matrix
             */
    //***********************************************************************
    static void InvertMatrix2(const MatrixType& InputMatrix,
                              MatrixType& InvertedMatrix,
                              TDataType& InputMatrixDet)
    {
        KRATOS_TRY
        InvertedMatrix.resize(2,2);

        InputMatrixDet = InputMatrix(0,0)*InputMatrix(1,1)-InputMatrix(0,1)*InputMatrix(1,0);

        InvertedMatrix(0,0) =  InputMatrix(1,1);
        InvertedMatrix(0,1) = -InputMatrix(0,1);
        InvertedMatrix(1,0) = -InputMatrix(1,0);
        InvertedMatrix(1,1) =  InputMatrix(0,0);

        InvertedMatrix/=InputMatrixDet;
        KRATOS_CATCH("")
    }

    //***********************************************************************
    /**calculates the determinant of a matrix of dimension 2*2
        (no check performed on matrix size)*/
    //***********************************************************************
    static inline TDataType Det2(const MatrixType& A)
    {
        return (A(0,0)*A(1,1)-A(0,1)*A(1,0));
    }

    //***********************************************************************
    /**calculates the determinant of a matrix of dimension 3*3
        (no check performed on matrix size)*/
    //***********************************************************************
    static inline TDataType Det3(const MatrixType& A)
    {
        //calculating the algebraic complements to the first line
        TDataType a = A(1,1)*A(2,2) - A(1,2)*A(2,1);
        TDataType b = A(1,0)*A(2,2) - A(1,2)*A(2,0);
        TDataType c = A(1,0)*A(2,1) - A(1,1)*A(2,0);

        return A(0,0)*a - A(0,1)*b + A(0,2)*c;
    }

    //***********************************************************************
    /**calculates the determinant of a matrix of dimension 2*2 or 3*3
        (no check performed on matrix size)*/
    //***********************************************************************
    static inline TDataType Det(const boost::numeric::ublas::bounded_matrix<TDataType,2,2>& A)
    {
        return (A(0,0)*A(1,1)-A(0,1)*A(1,0));
    }
    static inline TDataType Det(const boost::numeric::ublas::bounded_matrix<TDataType,3,3>& A)
    {
        //calculating the algebraic complements to the first line
        TDataType a = A(1,1)*A(2,2) - A(1,2)*A(2,1);
        TDataType b = A(1,0)*A(2,2) - A(1,2)*A(2,0);
        TDataType c = A(1,0)*A(2,1) - A(1,1)*A(2,0);

        return A(0,0)*a - A(0,1)*b + A(0,2)*c;
    }
    static inline TDataType Det(const MatrixType& A)
    {
        TDataType Det;
        if (A.size1()==1)
            Det = A(0,0);
        else if (A.size1()==2)
            Det = Det2(A);
        else
            Det = Det3(A);

        return Det;

    }

    //***********************************************************************
    /** performs the dot product of two vectors of dimension 3
        (no check performed on vector sizes)*/
    //***********************************************************************
    static inline TDataType Dot3(const VectorType& a, const VectorType& b)
    {
        return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    }

    //***********************************************************************
    /** performs the dot product of two vectors of arbitrary size
        (no check performed on vector sizes)*/
    //***********************************************************************
    static inline TDataType Dot(const VectorType& FirstVector, const VectorType& SecondVector)
    {
        typename VectorType::const_iterator i = FirstVector.begin();
        typename VectorType::const_iterator j = SecondVector.begin();
        TDataType temp = TDataType();
        while(i != FirstVector.end())
            temp += *i++ * *j++;
        return temp;
        //return std::inner_product(FirstVector.begin(), FirstVector.end(), SecondVector.begin(), TDataType());
    }

    //***********************************************************************
    /*
        calculates the norm of vector "a" which is assumed to be of size 3
        (no check is performed on the vector's size)
            */
    //***********************************************************************
    template<class TVectorType>
    static inline TDataType Norm3(const TVectorType& a)
    {
        TDataType temp = std::pow(a[0],2) + std::pow(a[1],2) + std::pow(a[2],2);
        temp = std::sqrt(temp);
        return temp;
    }

    //***********************************************************************
    /*
        calculates the norm of vector "a"
            */
    //***********************************************************************

    static inline TDataType Norm(const VectorType& a)
    {
        typename VectorType::const_iterator i = a.begin();
        TDataType temp = TDataType();
        while(i != a.end())
        {
            temp += (*i) * (*i);
            i++;
        }
        return sqrt(temp);
    }

    //***********************************************************************
    /*
        performs the vector product of the two input vectors a,b
        a,b are assumed to be of size 3 (no check is performed on vector sizes)
            */
    //***********************************************************************
    template<class T>
    static inline T CrossProduct(
        const T& a,
        const T& b
        )
    {
        T c(a);

        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];

        return c;
    }


    //identical but it assumes that the output vector is given already the correct size
    template< class T1, class T2 , class T3>
    static inline void CrossProduct(T1& c, const T2& a, const T3& b ){
        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];
    }

    static inline array_1d<TDataType, 3> UnitCrossProduct(const array_1d<TDataType, 3>& a, const array_1d<TDataType, 3>& b)
    {
        array_1d<TDataType, 3> c;

        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];

        const TDataType length = std::sqrt(inner_prod(c, c));
        c = (1.00/length) * c;
        return c;
    }

    /**
     * @brief Performs the unitary cross product of the two input vectors a,b
     * @details a,b are assumed to be of size 3 (no check is performed on vector sizes)
     * @param a First input vector
     * @param b Second input vector
     * @param c The resulting vector
     */
    template< class T1, class T2 , class T3>
    static inline void UnitCrossProduct(T1& c, const T2& a, const T3& b )
    {
        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];
        const double norm = norm_2(c);
        c/=norm;
    }

    //***********************************************************************
    /**
    returns a matrix :
        A = a.tensorproduct.b
        a,b are assumed to be of order 3, no check is performed on the size of the vectors
             */
    //***********************************************************************
    static inline MatrixType TensorProduct3(const VectorType& a, const VectorType& b)
    {
        MatrixType A(3,3);
        A(0,0)=a[0]*b[0];
        A(0,1)=a[0]*b[1];
        A(0,2)=a[0]*b[2];
        A(1,0)=a[1]*b[0];
        A(1,1)=a[1]*b[1];
        A(1,2)=a[1]*b[2];
        A(2,0)=a[2]*b[0];
        A(2,1)=a[2]*b[1];
        A(2,2)=a[2]*b[2];
        return A;
    }

    //***********************************************************************
    /**
        "InputMatrix" is ADDED to "Destination" matrix starting from
        InitialRow and InitialCol of the destination matrix
        "Destination" is assumed to be able to contain the "input matrix"
        (no check is performed on the bounds)
             */
    //***********************************************************************
    static inline void  AddMatrix(
        MatrixType& Destination,
        const MatrixType& InputMatrix,
        int InitialRow,
        int InitialCol)
    {
        KRATOS_TRY
        for(unsigned int i=0; i<InputMatrix.size1(); i++)
            for(unsigned int j=0; j<InputMatrix.size2(); j++)
                Destination(InitialRow+i, InitialCol+j) += InputMatrix(i,j);
        KRATOS_CATCH("")
    }
    //***********************************************************************
    /**
        "InputMatrix" is SUBTRACTED to "Destination" matrix starting from
        InitialRow and InitialCol of the destination matrix
        "Destination" is assumed to be able to contain the "input matrix"
        (no check is performed on the bounds)
             */
    //***********************************************************************
    static inline void  SubtractMatrix(
        MatrixType& Destination,
        const MatrixType& InputMatrix,
        int InitialRow,
        int InitialCol)
    {
        KRATOS_TRY
        for(unsigned int i=0; i<InputMatrix.size1(); i++)
            for(unsigned int j=0; j<InputMatrix.size2(); j++)
                Destination(InitialRow+i, InitialCol+j) -= InputMatrix(i,j);
        KRATOS_CATCH("")
    }

    //***********************************************************************
    /**
        "InputMatrix" is WRITTEN on "Destination" matrix starting from
        InitialRow and InitialCol of the destination matrix
        "Destination" is assumed to be able to contain the "input matrix"
        (no check is performed on the bounds)
        ATTENTION: Destination is overwritten!!
             */
    //***********************************************************************
    static inline void  WriteMatrix(
        MatrixType& Destination,
        const MatrixType& InputMatrix,
        int InitialRow,
        int InitialCol)
    {
        KRATOS_TRY
        for(unsigned int i=0; i<InputMatrix.size1(); i++)
            for(unsigned int j=0; j<InputMatrix.size2(); j++)
                Destination(InitialRow+i, InitialCol+j) = InputMatrix(i,j);
        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************
    //performs the Kroneker product of the Reduced MatrixType with the identity matrix of
    //size "dimension"
    static inline void  ExpandReducedMatrix(
        MatrixType& Destination,
        const MatrixType& ReducedMatrix,
        unsigned int dimension)
    {
        KRATOS_TRY
        unsigned int size=ReducedMatrix.size2();

        for (unsigned int i=0; i<size; i++)
        {
            unsigned int rowindex = i*dimension;
            for (unsigned int j=0; j<size; j++)
            {
                unsigned int colindex = j*dimension;
                for(unsigned int ii=0; ii<dimension; ii++)
                    Destination(rowindex+ii,colindex+ii) = ReducedMatrix(i,j);
            }
        }
        KRATOS_CATCH("")
    }
    //***********************************************************************
    //***********************************************************************
    //performs the Kroneker product of the Reduced MatrixType with the identity matrix of
    //size "dimension" ADDING to the destination matrix
    static inline void  ExpandAndAddReducedMatrix(
        MatrixType& Destination,
        const MatrixType& ReducedMatrix,
        unsigned int dimension)
    {
        KRATOS_TRY
        unsigned int size     = ReducedMatrix.size2();
        unsigned int rowindex = 0;
	unsigned int colindex = 0;
        for (unsigned int i=0; i<size; i++){
            rowindex = i*dimension;
            for (unsigned int j=0; j<size; j++){
                colindex = j*dimension;
                for(unsigned int ii=0; ii<dimension; ii++)
                    Destination(rowindex+ii,colindex+ii)+=ReducedMatrix(i,j);
            } }
        KRATOS_CATCH("")
    }

    //***********************************************************************
    /**
        performs x += coeff*y;
        no check on bounds is performed
             */
    //***********************************************************************
    static inline void  VecAdd(
        VectorType& x,
        const TDataType& coeff,
        const VectorType& y)
    {
        KRATOS_TRY
        unsigned int size=x.size();

        for (unsigned int i=0; i<size; i++)
        {
            x[i] += coeff * y[i];
        }
        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************

   /**
     * Transforms a stess vector into a matrix. Stresses are assumed to be stored
     * in the following way:
     * \f$ [ s11, s22, s33, s12, s23, s13 ] \f$ for 3D case and
     * \f$ [ s11, s22, s33, s12 ] \f$ for 2D case.
     * \f$ [ s11, s22, s12 ] \f$ for 2D case.
     * @param rStressVector the given stress vector
     * @return the corresponding stress tensor in matrix form
     */
    static inline MatrixType StressVectorToTensor(const VectorType& rStressVector)
    {
      KRATOS_TRY
      MatrixType StressTensor;

      if (rStressVector.size()==3)
        {
	  StressTensor.resize(2,2,false);
	  StressTensor(0,0) = rStressVector[0];
	  StressTensor(0,1) = rStressVector[2];
	  StressTensor(1,0) = rStressVector[2];
	  StressTensor(1,1) = rStressVector[1];
        }
      else if (rStressVector.size()==4)
        {
	  StressTensor.resize(3,3,false);
	  StressTensor(0,0) = rStressVector[0];
	  StressTensor(0,1) = rStressVector[3];
	  StressTensor(0,2) = 0.0;
	  StressTensor(1,0) = rStressVector[3];
	  StressTensor(1,1) = rStressVector[1];
	  StressTensor(1,2) = 0.0;
	  StressTensor(2,0) = 0.0;
	  StressTensor(2,1) = 0.0;
	  StressTensor(2,2) = rStressVector[2];
	}
      else if (rStressVector.size()==6)
        {
	  StressTensor.resize(3,3,false);
	  StressTensor(0,0) = rStressVector[0];
	  StressTensor(0,1) = rStressVector[3];
	  StressTensor(0,2) = rStressVector[5];
	  StressTensor(1,0) = rStressVector[3];
	  StressTensor(1,1) = rStressVector[1];
	  StressTensor(1,2) = rStressVector[4];
	  StressTensor(2,0) = rStressVector[5];
	  StressTensor(2,1) = rStressVector[4];
	  StressTensor(2,2) = rStressVector[2];
        }

      return StressTensor;
      KRATOS_CATCH("")

    }
    

    //***********************************************************************
    //***********************************************************************

   /**
     * Transforms a  vector into a symmetric matrix.  Components are assumed to be stored
     * in the following way:
     * \f$ [ s11, s22, s33, s12, s23, s13 ] \f$ for 3D case and
     * \f$ [ s11, s22, s33, s12 ] \f$ for 2D case.
     * \f$ [ s11, s22, s12 ] \f$ for 2D case.
     * @param rVector the given stress vector
     * @return the corresponding Tensor in matrix form
     */
    static inline MatrixType VectorToSymmetricTensor(const VectorType& rVector)
    {
      KRATOS_TRY
      MatrixType Tensor;

      if (rVector.size()==3)
        {
	  Tensor.resize(2,2,false);
	  Tensor(0,0) = rVector[0];
	  Tensor(0,1) = rVector[2];
	  Tensor(1,0) = rVector[2];
	  Tensor(1,1) = rVector[1];
        }
      else if (rVector.size()==4)
        {
	  Tensor.resize(3,3,false);
	  Tensor(0,0) = rVector[0];
	  Tensor(0,1) = rVector[3];
	  Tensor(0,2) = 0.0;
	  Tensor(1,0) = rVector[3];
	  Tensor(1,1) = rVector[1];
	  Tensor(1,2) = 0.0;
	  Tensor(2,0) = 0.0;
	  Tensor(2,1) = 0.0;
	  Tensor(2,2) = rVector[2];
	}
      else if (rVector.size()==6)
        {
	  Tensor.resize(3,3,false);
	  Tensor(0,0) = rVector[0];
	  Tensor(0,1) = rVector[3];
	  Tensor(0,2) = rVector[5];
	  Tensor(1,0) = rVector[3];
	  Tensor(1,1) = rVector[1];
	  Tensor(1,2) = rVector[4];
	  Tensor(2,0) = rVector[5];
	  Tensor(2,1) = rVector[4];
	  Tensor(2,2) = rVector[2];
        }

      return Tensor;
      KRATOS_CATCH("")

    }
    
    
    //***********************************************************************
    //***********************************************************************
    /// sign function
    static inline int Sign(const TDataType& ThisDataType)
    {
        KRATOS_TRY
        const TDataType& x = ThisDataType;
        return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
        KRATOS_CATCH("")
    }
    

    /**
     * Transforms a strain vector into a matrix. Strains are assumed to be stored
     * in the following way:
     * \f$ [ e11, e22, e33, 2*e12, 2*e23, 2*e13 ] \f$ for 3D case and
     * \f$ [ e11, e22, e33, 2*e12 ] \f$ for 2D case.
     * \f$ [ e11, e22, 2*e12 ] \f$ for 2D case.
     * Hence the deviatoric components of the strain vector are divided by 2
     * while they are stored into the matrix
     * @param rStrainVector the given strain vector
     * @return the corresponding strain tensor in matrix form
     */
    static inline MatrixType StrainVectorToTensor( const VectorType& rStrainVector)
    {
      KRATOS_TRY
      MatrixType StrainTensor;

      if (rStrainVector.size()==3)
        {
	  StrainTensor.resize(2,2, false);
	  StrainTensor(0,0) = rStrainVector[0];
	  StrainTensor(0,1) = 0.5*rStrainVector[2];
	  StrainTensor(1,0) = 0.5*rStrainVector[2];
	  StrainTensor(1,1) = rStrainVector[1];
        }
      else if (rStrainVector.size()==4)
        {
	  StrainTensor.resize(3,3, false);
	  StrainTensor(0,0) = rStrainVector[0];
	  StrainTensor(0,1) = 0.5*rStrainVector[3];
	  StrainTensor(0,2) = 0;
	  StrainTensor(1,0) = 0.5*rStrainVector[3];
	  StrainTensor(1,1) = rStrainVector[1];
	  StrainTensor(1,2) = 0;
	  StrainTensor(2,0) = 0;
	  StrainTensor(2,1) = 0;
	  StrainTensor(2,2) = rStrainVector[2];
        }
      else if (rStrainVector.size()==6)
        {
	  StrainTensor.resize(3,3, false);
	  StrainTensor(0,0) = rStrainVector[0];
	  StrainTensor(0,1) = 0.5*rStrainVector[3];
	  StrainTensor(0,2) = 0.5*rStrainVector[5];
	  StrainTensor(1,0) = 0.5*rStrainVector[3];
	  StrainTensor(1,1) = rStrainVector[1];
	  StrainTensor(1,2) = 0.5*rStrainVector[4];
	  StrainTensor(2,0) = 0.5*rStrainVector[5];
	  StrainTensor(2,1) = 0.5*rStrainVector[4];
	  StrainTensor(2,2) = rStrainVector[2];
        }

      return StrainTensor;
      KRATOS_CATCH("")
    }


    /**
    * Transforms a given symmetric Strain Tensor to Voigt Notation:
    * in the 3D case: from a second order tensor (3*3) MatrixType  to a corresponing (6*1) VectorType
    * \f$ [ e11, e22, e33, 2*e12, 2*e23, 2*e13 ] \f$ for 3D case and
    * in the 2D case: from a second order tensor (3*3) MatrixType  to a corresponing (4*1) VectorType
    * \f$ [ e11, e22, e33, 2*e12 ] \f$ fir 2D case.
    * in the 2D case: from a second order tensor (2*2) MatrixType  to a corresponing (3*1) VectorType
    * \f$ [ e11, e22, 2*e12 ] \f$ fir 2D case.
    * @param rStrainTensor the given symmetric second order strain tensor
    * @return the corresponding strain tensor in vector form
    */

    static inline VectorType StrainTensorToVector( const MatrixType& rStrainTensor, unsigned int rSize = 0 )
    {
      KRATOS_TRY

      VectorType StrainVector;

     if(rSize == 0){
	if(rStrainTensor.size1() == 2)
	  rSize = 3;
	if(rStrainTensor.size1() == 3)
	  rSize = 6;
      }

      if (rSize == 3)
        {
	  StrainVector.resize(3);
	  StrainVector[0] = rStrainTensor(0,0);
	  StrainVector[1] = rStrainTensor(1,1);
	  StrainVector[2] = 2.00*rStrainTensor(0,1);
        }
      else if (rSize == 4)
        {
	  StrainVector.resize(4);
	  StrainVector[0] = rStrainTensor(0,0);
	  StrainVector[1] = rStrainTensor(1,1);
	  StrainVector[2] = rStrainTensor(2,2);
	  StrainVector[3] = 2.00*rStrainTensor(0,1);
        }
      else if (rSize == 6)
        {
	  StrainVector.resize(6);
	  StrainVector[0] = rStrainTensor(0,0);
	  StrainVector[1] = rStrainTensor(1,1);
	  StrainVector[2] = rStrainTensor(2,2);
	  StrainVector[3] = 2.00*rStrainTensor(0,1);
	  StrainVector[4] = 2.00*rStrainTensor(1,2);
	  StrainVector[5] = 2.00*rStrainTensor(0,2);
        }


      return StrainVector;
      KRATOS_CATCH("")
     }


    /**
    * Transforms a given symmetric Stress Tensor to Voigt Notation:
    * in the 3D case: from a second order tensor (3*3) MatrixType  to a corresponing (6*1) VectorType
    * in the 3D case: from a second order tensor (3*3) MatrixType  to a corresponing (4*1) VectorType
    * in the 2D case: from a second order tensor (2*2) MatrixType  to a corresponing (3*1) VectorType
    * @param rStressTensor the given symmetric second order stress tensor
    * @return the corresponding stress tensor in vector form
    */
    static inline VectorType StressTensorToVector(const MatrixType& rStressTensor, unsigned int rSize = 0)
    {

      KRATOS_TRY

      VectorType StressVector;
      
      if(rSize == 0){
	if(rStressTensor.size1() == 2)
	  rSize = 3;
	if(rStressTensor.size1() == 3)
	  rSize = 6;
      }

      if (rSize==3)
        {
	  StressVector.resize(3);
	  StressVector[0]= rStressTensor(0,0);
	  StressVector[1]= rStressTensor(1,1);
	  StressVector[2]= rStressTensor(0,1);
        }
      else if (rSize==4)
        {
	  StressVector.resize(4);
	  StressVector[0]= rStressTensor(0,0);
	  StressVector[1]= rStressTensor(1,1);
	  StressVector[2]= rStressTensor(2,2);
	  StressVector[3]= rStressTensor(0,1);
        }
      else if (rSize==6)
        {
	  StressVector.resize(6);
	  StressVector[0]= rStressTensor(0,0);
	  StressVector[1]= rStressTensor(1,1);
	  StressVector[2]= rStressTensor(2,2);
	  StressVector[3]= rStressTensor(0,1);
	  StressVector[4]= rStressTensor(1,2);
	  StressVector[5]= rStressTensor(0,2);
        }

        
      return StressVector;
      KRATOS_CATCH("")
     }



    /**
    * Transforms a given symmetric Tensor to Voigt Notation:
    * in the 3D case: from a second order tensor (3*3) MatrixType  to a corresponing (6*1) VectorType
    * in the 3D case: from a second order tensor (3*3) MatrixType  to a corresponing (4*1) VectorType
    * in the 2D case: from a second order tensor (2*2) MatrixType  to a corresponing (3*1) VectorType
    * @param rStressTensor the given symmetric second order stress tensor
    * @return the corresponding stress tensor in vector form
    */
    static inline VectorType SymmetricTensorToVector(const MatrixType& rTensor, unsigned int rSize = 0)
    {

      KRATOS_TRY

      VectorType vector;
      
      if(rSize == 0){
	if(rTensor.size1() == 2)
	  rSize = 3;
	if(rTensor.size1() == 3)
	  rSize = 6;
      }

      if (rSize==3)
        {
	  vector.resize(3);
	  vector[0]= rTensor(0,0);
	  vector[1]= rTensor(1,1);
	  vector[2]= rTensor(0,1);
        }
      else if (rSize==4)
        {
	  vector.resize(4);
	  vector[0]= rTensor(0,0);
	  vector[1]= rTensor(1,1);
	  vector[2]= rTensor(2,2);
	  vector[3]= rTensor(0,1);
        }
      else if (rSize==6)
        {
	  vector.resize(6);
	  vector[0]= rTensor(0,0);
	  vector[1]= rTensor(1,1);
	  vector[2]= rTensor(2,2);
	  vector[3]= rTensor(0,1);
	  vector[4]= rTensor(1,2);
	  vector[5]= rTensor(0,2);
        }

        
      return vector;
      KRATOS_CATCH("")
     }

    /**
     * Solve a1 * x + b1 * y = c1
     *       a2 * x + b2 * y = c2
     */
    static void Solve(
        const TDataType& a1, const TDataType& b1, const TDataType& c1,
        const TDataType& a2, const TDataType& b2, const TDataType& c2,
        TDataType& x, TDataType& y
    )
    {
        x = (c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1);
        y = (c2 * a1 - c1 * a2) / (a1 * b2 - a2 * b1);
    }

    /*@} */
    /**@name Acces */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */


    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Acces */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    MathUtils(void);

    MathUtils(MathUtils& rSource);


    /*@} */

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_MATH_UTILS  defined */

