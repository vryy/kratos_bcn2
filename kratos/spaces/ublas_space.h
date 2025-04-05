//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#if !defined(KRATOS_UBLAS_SPACE_H_INCLUDED )
#define  KRATOS_UBLAS_SPACE_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <iomanip>
#include <cstddef>
#include <numeric>
#include <type_traits>

#ifdef _OPENMP
#include "omp.h"
#endif

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "includes/matrix_market_interface.h"
#include "utilities/openmp_utils.h"


namespace Kratos
{
#ifdef _OPENMP
// The function object multiplies an element by a Factor

template <class Type>
class MultValueNoAdd
{
private:
    Type Factor; // The value to multiply by
public:
    // Constructor initializes the value to multiply by

    MultValueNoAdd(const Type& _Val) : Factor(_Val)
    {
    }

    // The function call for the element to be multiplied

    inline Type operator () (const Type& elem) const
    {
        return elem * Factor;
    }
};

template <class Type>
class MultAndAddValue
{
private:
    Type Factor; // The value to multiply by
public:
    // Constructor initializes the value to multiply by

    MultAndAddValue(const Type& _Val) : Factor(_Val)
    {
    }

    // The function call for the element to be multiplied

    inline Type operator () (const Type& elem1, const Type& elem2) const
    {
        return elem1 * Factor + elem2;
    }
};
#endif

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
 * @class UblasSpace
 * @ingroup KratosCore
 * @brief A class template for handling data types, matrices, and vectors in a Ublas space.
 * @details This class template is designed to work with different data types, matrix types, and vector types
 * within a Ublas space. It provides typedefs and utilities for managing these types effectively.
 * @tparam TDataType The data type used in the Ublas space.
 * @tparam TMatrixType The matrix type used in the Ublas space.
 * @tparam TVectorType The vector type used in the Ublas space.
 * @author Riccardo Rossi
 */
template<class TDataType, class TMatrixType, class TVectorType>
class UblasSpace
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of UblasSpace
    KRATOS_CLASS_POINTER_DEFINITION(UblasSpace);

    typedef TDataType DataType;

    typedef TMatrixType MatrixType;

    typedef TVectorType VectorType;

    typedef KRATOS_INDEX_TYPE IndexType;

    typedef KRATOS_SIZE_TYPE SizeType;

    typedef typename boost::shared_ptr< TMatrixType > MatrixPointerType;
    typedef typename boost::shared_ptr< TVectorType > VectorPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UblasSpace()
    {
    }

    /// Destructor.
    virtual ~UblasSpace()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static MatrixPointerType CreateEmptyMatrixPointer()
    {
        return MatrixPointerType(new TMatrixType());
    }

    static VectorPointerType CreateEmptyVectorPointer()
    {
        return VectorPointerType(new TVectorType());
    }

    /// return size of vector rV
    static IndexType Size(VectorType const& rV)
    {
        return rV.size();
    }

    /// return number of rows of rM
    static IndexType Size1(MatrixType const& rM)
    {
        return rM.size1();
    }

    /// return number of columns of rM
    static IndexType Size2(MatrixType const& rM)
    {
        return rM.size2();
    }

    /// rXi = rMij
    static void GetColumn(unsigned int j, const MatrixType& rM, VectorType& rX)
    {
        noalias(rX) = column(rM, j);
    }

    ///////////////////////////////// TODO: Take a close look to this method!!!!!!!!!!!!!!!!!!!!!!!!!
    /// rMij = rXi
    static void SetColumn(unsigned int j, MatrixType& rM, const VectorType& rX)
    {
        noalias(column(rM, j)) = rX;
    }

    /// rY = rX
    static void Copy(MatrixType const& rX, MatrixType& rY)
    {
        rY.assign(rX);
    }

    /// rY = rX^T
    static void TransposeCopy(MatrixType const& rX, MatrixType& rY)
    {
        rY.assign(trans(rX));
    }

    /// rY = rX
    static void Copy(VectorType const& rX, VectorType& rY)
    {
#ifndef _OPENMP
        rY.assign(rX);
#else

        const int size = rX.size();
        if (rY.size() != static_cast<unsigned int>(size))
            rY.resize(size, false);

        #pragma omp parallel for
        for (int i = 0; i < size; i++)
            rY[i] = rX[i];
#endif
    }

    /// rX * rY
    static TDataType Dot(VectorType const& rX, VectorType const& rY)
    {
#ifndef _OPENMP
        return inner_prod(rX, rY);
#else
        std::vector<unsigned int> partition;
        int number_of_threads = omp_get_max_threads();
        OpenMPUtils::CreatePartition(number_of_threads, rX.size(), partition);

        vector< TDataType > partial_results(number_of_threads);

        #pragma omp parallel for
        for (int i = 0; i < number_of_threads; i++)
        {
            partial_results[i] = std::inner_product(rX.data().begin() + partition[i],
                                                    rX.data().begin() + partition[i + 1],
                                                    rY.data().begin() + partition[i],
                                                    TDataType());
        }

        auto total = TDataType();
        for (int i = 0; i < number_of_threads; i++)
            total += partial_results[i];
        return total;
#endif
    }

    /// ||rX||2
    static TDataType TwoNorm(VectorType const& rX)
    {
        return sqrt(Dot(rX, rX));
    }

    // Frobenious norm
    static TDataType TwoNorm(const MatrixType& rA)
    {
        if constexpr (std::is_same<MatrixType, Matrix>::value)
        {
            TDataType aux_sum = TDataType();
            #pragma omp parallel for reduction(+:aux_sum)
            for (int i=0; i<static_cast<int>(rA.size1()); ++i) {
                for (int j=0; j<static_cast<int>(rA.size2()); ++j) {
                    aux_sum += std::pow(rA(i,j),2);
                }
            }
            return std::sqrt(aux_sum);
        }
        else if constexpr (std::is_same<MatrixType, compressed_matrix<TDataType> >::value)
        {
            TDataType aux_sum = TDataType();

            const auto& r_values = rA.value_data();

            #pragma omp parallel for reduction(+:aux_sum)
            for (int i=0; i<static_cast<int>(r_values.size()); ++i) {
                aux_sum += std::pow(r_values[i] , 2);
            }
            return std::sqrt(aux_sum);
        }
        else
            KRATOS_ERROR << __FUNCTION__ << " is not implemented";
    }

    /**
     * This method computes the Jacobi norm
     * @param rA The matrix to compute the Jacobi norm
     * @return aux_sum: The Jacobi norm
     */
    static TDataType JacobiNorm(const MatrixType& rA)
    {
        if constexpr (std::is_same<MatrixType, Matrix>::value)
        {
            TDataType aux_sum = TDataType();
            #pragma omp parallel for reduction(+:aux_sum)
            for (int i=0; i<static_cast<int>(rA.size1()); ++i) {
                for (int j=0; j<static_cast<int>(rA.size2()); ++j) {
                    if (i != j) {
                        aux_sum += std::abs(rA(i,j));
                    }
                }
            }
            return aux_sum;
        }
        else if constexpr (std::is_same<MatrixType, compressed_matrix<TDataType> >::value)
        {
            TDataType aux_sum = TDataType();

            typedef typename MatrixType::const_iterator1 t_it_1;
            typedef typename MatrixType::const_iterator2 t_it_2;

            for (t_it_1 it_1 = rA.begin1(); it_1 != rA.end1(); ++it_1) {
                for (t_it_2 it_2 = it_1.begin(); it_2 != it_1.end(); ++it_2) {
                    if (it_2.index1() != it_2.index2()) {
                        aux_sum += std::abs(*it_2);
                    }
                }
            }

            return aux_sum;
        }
        else
            KRATOS_ERROR << __FUNCTION__ << " is not implemented";
    }

    static void Mult(Matrix& rA, VectorType& rX, VectorType& rY)
    {
        axpy_prod(rA, rX, rY, true);
    }

    static void Mult(compressed_matrix<TDataType>& rA, VectorType& rX, VectorType& rY)
    {
#ifndef _OPENMP
        axpy_prod(rA, rX, rY, true);
#else
        ParallelProductNoAdd(rA, rX, rY);
#endif
    }

    static void TransposeMult(MatrixType& rA, VectorType& rX, VectorType& rY)
    {
        axpy_prod(rX, rA, rY, true);
    } // rY = rAT * rX

    static inline SizeType GraphDegree(IndexType i, TMatrixType& A)
    {
        typename MatrixType::iterator1 a_iterator = A.begin1();
        std::advance(a_iterator, i);
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
        return ( std::distance(a_iterator.begin(), a_iterator.end()));
#else
        return ( std::distance(begin(a_iterator, boost::numeric::ublas::iterator1_tag()),
                               end(a_iterator, boost::numeric::ublas::iterator1_tag())));
#endif
    }

    static inline void GraphNeighbors(IndexType i, TMatrixType& A, std::vector<IndexType>& neighbors)
    {
        neighbors.clear();
        typename MatrixType::iterator1 a_iterator = A.begin1();
        std::advance(a_iterator, i);
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
        for (typename MatrixType::iterator2 row_iterator = a_iterator.begin();
                row_iterator != a_iterator.end(); ++row_iterator)
        {
#else
        for (typename MatrixType::iterator2 row_iterator = begin(a_iterator,
                boost::numeric::ublas::iterator1_tag());
                row_iterator != end(a_iterator,
                                    boost::numeric::ublas::iterator1_tag()); ++row_iterator)
        {
#endif
            neighbors.push_back(row_iterator.index2());
        }
    }

    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise
    static void InplaceMult(VectorType& rX, const TDataType A)
    {

        if (A == 1.00)
        {
        }
        else if (A == -1.00)
        {
#ifndef _OPENMP
            typename VectorType::iterator x_iterator = rX.begin();
            typename VectorType::iterator end_iterator = rX.end();
            while (x_iterator != end_iterator)
            {
                *x_iterator = -*x_iterator;
                x_iterator++;

            }
#else
            const int size = rX.size();

            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] = -rX[i];

#endif
        }
        else
        {
#ifndef _OPENMP
            rX *= A;
#else
            const int size = rX.size();

            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] *= A;
#endif
        }
    }

    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise
    //ATTENTION it is assumed no aliasing between rX and rY
    // X = A*y;
    static void Assign(VectorType& rX, const TDataType A, const VectorType& rY)
    {
#ifndef _OPENMP
        if (A == 1.00)
            noalias(rX) = rY;
        else if (A == -1.00)
            noalias(rX) = -rY;
        else
            noalias(rX) = A*rY;
#else
        const int size = rY.size();
        if (rX.size() != static_cast<unsigned int>(size) )
            rX.resize(size, false);

        if (A == 1.00)
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] = rY[i];
        }
        else if (A == -1.00)
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] = -rY[i];
        }
        else
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] = A * rY[i];
        }
#endif
    }

    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise
    //ATTENTION it is assumed no aliasing between rX and rY
    // X += A*y;
    static void UnaliasedAdd(VectorType& rX, const TDataType A, const VectorType& rY)
    {
#ifndef _OPENMP
        if (A == 1.00)
            noalias(rX) += rY;
        else if (A == -1.00)
            noalias(rX) -= rY;
        else
            noalias(rX) += A*rY;
#else
        const int size = rY.size();
        if (rX.size() != static_cast<unsigned int>(size) )
            rX.resize(size, false);

        if (A == 1.00)
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] += rY[i];
        }
        else if (A == -1.00)
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] -= rY[i];
        }
        else
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] += A * rY[i];
        }
#endif
    }

    //********************************************************************

    static void ScaleAndAdd(const TDataType A, const VectorType& rX, const TDataType B, const VectorType& rY, VectorType& rZ) // rZ = (A * rX) + (B * rY)
    {
        Assign(rZ, A, rX); //rZ = A*rX
        UnaliasedAdd(rZ, B, rY); //rZ += B*rY
    }

    static void ScaleAndAdd(const TDataType A, const VectorType& rX, const TDataType B, VectorType& rY) // rY = (A * rX) + (B * rY)
    {
        InplaceMult(rY, B);
        UnaliasedAdd(rY, A, rX);
    }


    /// rA[i] * rX
    static TDataType RowDot(unsigned int i, MatrixType& rA, VectorType& rX)
    {
        return inner_prod(row(rA, i), rX);
    }

    static void SetValue(VectorType& rX, IndexType i, TDataType value)
    {
        rX[i] = value;
    }

    /// rX = A
    static void Set(VectorType& rX, TDataType A)
    {
        std::fill(rX.begin(), rX.end(), A);
    }

    static void Resize(MatrixType& rA, SizeType m, SizeType n)
    {
        rA.resize(m, n, false);
    }

    static void Resize(VectorType& rX, SizeType n)
    {
        rX.resize(n, false);
    }

    static void Clear(MatrixPointerType& pA)
    {
        pA->clear();
        pA->resize(0, 0, false);
    }

    static void Clear(VectorPointerType& pX)
    {
        pX->clear();
        pX->resize(0, false);
    }

    template<class TOtherMatrixType>
    inline static void ClearData(TOtherMatrixType& rA)
    {
        rA.clear();
    }

    inline static void ClearData(compressed_matrix<TDataType>& rA)
    {
        rA.clear();
    }

    inline static void ClearData(VectorType& rX)
    {
        rX = VectorType();
    }

    template<class TOtherMatrixType>
    inline static void ResizeData(TOtherMatrixType& rA, SizeType m)
    {
        rA.resize(m, false);
#ifndef _OPENMP
        std::fill(rA.begin(), rA.end(), TDataType());
#else
        ParallelFill(rA.begin(), rA.end(), TDataType());

#endif
    }

    inline static void ResizeData(compressed_matrix<TDataType>& rA, SizeType m)
    {
        rA.value_data().resize(m);
#ifndef _OPENMP
        std::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());
#else
        ParallelFill(rA.value_data().begin(), rA.value_data().end(), TDataType());
#endif
    }

    inline static void ResizeData(VectorType& rX, SizeType m)
    {
        rX.resize(m, false);
#ifndef _OPENMP
        std::fill(rX.begin(), rX.end(), TDataType());
#else
        ParallelFill(rX.begin(), rX.end(), TDataType());
#endif
    }

    template<class TOtherMatrixType>
    inline static void SetToZero(TOtherMatrixType& rA)
    {
#ifndef _OPENMP
        std::fill(rA.begin(), rA.end(), TDataType());
#else
        ParallelFill(rA.begin(), rA.end(), TDataType());
#endif
    }

    inline static void SetToZero(compressed_matrix<TDataType>& rA)
    {
#ifndef _OPENMP
        std::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());
#else
        ParallelFill(rA.value_data().begin(), rA.value_data().end(), TDataType());
#endif
    }

    inline static void SetToZero(VectorType& rX)
    {
#ifndef _OPENMP
        std::fill(rX.begin(), rX.end(), TDataType());
#else
        ParallelFill(rX.begin(), rX.end(), TDataType());
#endif
    }

    template<class TOtherMatrixType, class TEquationIdVectorType>
    inline static void AssembleLHS(
        MatrixType& A,
        TOtherMatrixType& LHS_Contribution,
        TEquationIdVectorType& EquationId
    )
    {
        unsigned int system_size = A.size1();
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];
            if (i_global < system_size)
            {
                for (unsigned int j_local = 0; j_local < local_size; j_local++)
                {
                    unsigned int j_global = EquationId[j_local];
                    if (j_global < system_size)
                        A(i_global, j_global) += LHS_Contribution(i_local, j_local);
                }
            }
        }
    }

    ///@}
    ///@name Access
    ///@{

    //***********************************************************************

    inline static TDataType GetValue(const VectorType& x, std::size_t I)
    {
        return x[I];
    }

    //***********************************************************************

    static void GatherValues(const VectorType& x, const std::vector<std::size_t>& IndexArray, TDataType* pValues)
    {
        KRATOS_TRY

        for (std::size_t i = 0; i < IndexArray.size(); i++)
            pValues[i] = x[IndexArray[i]];

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Inquiry
    ///@{

    inline static constexpr bool IsDistributed()
    {
        return false;
    }

    /// Print the information about the matrix. Depending on the level, different information will be printed
    static void PrintMatrixInfo(std::ostream& rOStream, const MatrixType& rA, const int level)
    {
        typedef typename MatrixType::const_iterator1 t_it_1;
        typedef typename MatrixType::const_iterator2 t_it_2;

        if (level > 0)
        {
            rOStream << "Number of rows: " << Size1(rA)
                     << ", columns: " << Size2(rA);
        }

        if (level > 1)
        {
            std::size_t num_entries = 0;

            if constexpr (std::is_same<MatrixType, Matrix>::value)
            {
                for (IndexType i = 0; i < Size1(rA); ++i)
                    for (IndexType j = 0; j < Size2(rA); ++j)
                        if (rA(i, j) != 0) ++num_entries;
            }
            else if constexpr (std::is_same<MatrixType, compressed_matrix<TDataType> >::value)
            {
                for (t_it_1 it_1 = rA.begin1(); it_1 != rA.end1(); ++it_1)
                    for (t_it_2 it_2 = it_1.begin(); it_2 != it_1.end(); ++it_2)
                        if (*it_2 != 0.0) ++num_entries;
            }

            rOStream << ", entries: " << num_entries;
        }

        if (level > 2)
        {
            DataType max_diag = -1e99, sum_diag = 0.0;

            if constexpr (std::is_same<MatrixType, Matrix>::value)
            {
                for (IndexType i = 0; i < Size1(rA); ++i)
                    for (IndexType j = 0; j < Size2(rA); ++j)
                        if (i == j)
                        {
                            const DataType v = std::abs(rA(i, j));
                            sum_diag += v;
                            if (v > max_diag)
                                max_diag = v;
                        }
            }
            else if constexpr (std::is_same<MatrixType, compressed_matrix<TDataType> >::value)
            {
                for (t_it_1 it_1 = rA.begin1(); it_1 != rA.end1(); ++it_1)
                    for (t_it_2 it_2 = it_1.begin(); it_2 != it_1.end(); ++it_2)
                        if (it_2.index1() == it_2.index2())
                        {
                            const DataType v = std::abs(*it_2);
                            sum_diag += v;
                            if (v > max_diag)
                                max_diag = v;
                        }
            }

            rOStream << std::scientific << std::setprecision(6)
                     << ", max of abs of diag: " << max_diag
                     << ", sum of abs of diag: " << sum_diag;
        }

        rOStream << std::endl;
    }

    /// Print the information about the vector. Depending on the level, different information will be printed
    static void PrintVectorInfo(std::ostream& rOStream, const VectorType& rX, const int level)
    {
        typedef typename VectorType::const_iterator t_it;

        if (level > 0)
        {
            rOStream << "Size: " << Size(rX);
        }

        if (level > 1)
        {
            std::size_t num_entries = 0;
            for (t_it it = rX.begin(); it != rX.end(); ++it)
            {
                if (*it != 0.0)
                    ++num_entries;
            }
            rOStream << ", entries: " << num_entries;
        }

        if (level > 2)
        {
            DataType max_entry = -1e99, sum_entry = 0.0;
            for (t_it it = rX.begin(); it != rX.end(); ++it)
            {
                const DataType v = std::abs(*it);
                sum_entry += v;
                if (v > max_entry)
                    max_entry = v;
            }
            rOStream << std::scientific << std::setprecision(6)
                     << ", max of abs of entry: " << max_entry
                     << ", sum of abs of entry: " << sum_entry;
        }

        rOStream << std::endl;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "UBlasSpace";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "UBlasSpace";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    //***********************************************************************

    template< class TOtherMatrixType >
    static bool WriteMatrixMarketMatrix(const char *FileName, TOtherMatrixType &M, bool Symmetric)
    {
        // Use full namespace in call to make sure we are not calling this function recursively
        return Kratos::WriteMatrixMarketMatrix(FileName,M,Symmetric);
    }

    template< class VectorType >
    static bool WriteMatrixMarketVector(const char *FileName, VectorType& V)
    {
        // Use full namespace in call to make sure we are not calling this function recursively
        return Kratos::WriteMatrixMarketVector(FileName,V);
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

#ifdef _OPENMP
    //y += A*x in parallel
    static void ParallelProductNoAdd(const MatrixType& A, const VectorType& in, VectorType& out)
    {
        //create partition
        std::vector<unsigned int> partition;
        unsigned int number_of_threads = omp_get_max_threads();
        unsigned int number_of_initialized_rows = A.filled1() - 1;
        OpenMPUtils::CreatePartition(number_of_threads, number_of_initialized_rows, partition);
        //parallel loop
        #pragma omp parallel
        {
            int thread_id = omp_get_thread_num();
            int number_of_rows = partition[thread_id + 1] - partition[thread_id];
            typename compressed_matrix<TDataType>::index_array_type::const_iterator row_iter_begin = A.index1_data().begin() + partition[thread_id];
            typename compressed_matrix<TDataType>::index_array_type::const_iterator index_2_begin = A.index2_data().begin()+*row_iter_begin;
            typename compressed_matrix<TDataType>::value_array_type::const_iterator value_begin = A.value_data().begin()+*row_iter_begin;
            //                  typename VectorType::iterator output_vec_begin = out.begin()+partition[thread_id];


            partial_product_no_add(number_of_rows,
                                   row_iter_begin,
                                   index_2_begin,
                                   value_begin,
                                   in,
                                   partition[thread_id],
                                   out
                                  );
        }
    }

    /**
     * calculates partial product resetting to Zero the output before
     */
    static void partial_product_no_add(
        int number_of_rows,
        typename compressed_matrix<TDataType>::index_array_type::const_iterator row_begin,
        typename compressed_matrix<TDataType>::index_array_type::const_iterator index2_begin,
        typename compressed_matrix<TDataType>::value_array_type::const_iterator value_begin,
        const VectorType& input_vec,
        unsigned int output_begin_index,
        VectorType& output_vec
        //                 typename VectorType::iterator output_vec_begin
    )
    {
        int row_size;
        int kkk = output_begin_index;
        typename MatrixType::index_array_type::const_iterator row_it = row_begin;
        for (int k = 0; k < number_of_rows; k++)
        {
            row_size = *(row_it + 1)-*row_it;
            row_it++;
            TDataType t = TDataType();

            for (int i = 0; i < row_size; i++)
                t += *value_begin++ * (input_vec[*index2_begin++]);

            output_vec[kkk++] = t;
            //                 *output_vec_begin++ = t;

        }
    }
#endif

    template <class TIterartorType>
    static void ParallelFill(TIterartorType Begin, TIterartorType End, TDataType const& Value)
    {
#ifndef _OPENMP
        std::fill(Begin, End, Value);
#else
        int size = End-Begin;
        //#pragma omp parallel for
        for (int i = 0; i < size; i++)
        {
            *(Begin+i) = Value;
        }
#endif
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

    /// Assignment operator.
    UblasSpace & operator=(UblasSpace const& rOther);

    /// Copy constructor.
    UblasSpace(UblasSpace const& rOther);

    ///@}
}; // Class UblasSpace

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template<class TDataType, class TMatrixType, class TVectorType>
inline std::ostream& operator << (std::ostream& rOStream,
             const UblasSpace<TDataType, TMatrixType, TVectorType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_UBLAS_SPACE_H_INCLUDED  defined
