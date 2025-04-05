//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

#if !defined(KRATOS_NODE_H_INCLUDED )
#define  KRATOS_NODE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes

// Project includes
#include "includes/define.h"
#include "geometries/point.h"
#include "includes/dof.h"
#include "containers/pointer_vector_set.h"
#include "containers/data_value_container.h"
#include "containers/variables_list_data_value_container.h"
#include "utilities/indexed_object.h"
#include "containers/flags.h"

#include "containers/weak_pointer_vector.h"

#ifdef _OPENMP
#include "omp.h"
#endif


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

/// This class defines the node
/** The node class from Kratos is defined in this class
*/
template<std::size_t TDimension, class TDofType = Dof<KRATOS_DOUBLE_TYPE> >
class Node : public Point<TDimension, typename TDofType::DataType>,  public IndexedObject, public Flags
{
    class GetDofKey
    {
    public:

        using argument_type = TDofType;
        using result_type = VariableData::KeyType;

        VariableData::KeyType operator()(TDofType const & This) const
        {
            return This.GetVariable().Key();
        }
    };

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Node
    KRATOS_CLASS_POINTER_DEFINITION(Node);

    typedef typename TDofType::DataType DataType;

    typedef Node<TDimension, TDofType> NodeType;

    typedef Point<TDimension, DataType> BaseType;

    typedef Point<TDimension, DataType> PointType;

    typedef TDofType DofType;

    typedef KRATOS_INDEX_TYPE IndexType;

    typedef KRATOS_SIZE_TYPE SizeType;

    typedef PointerVectorSet<TDofType, GetDofKey> DofsContainerType;

    typedef VariablesListDataValueContainer SolutionStepsNodalDataContainerType;

    typedef VariablesListDataValueContainer::BlockType BlockType;

    typedef Variable<DataType> DoubleVariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Node()
        : BaseType()
        , IndexedObject(0)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition()
    {
        CreateSolutionStepData();

#ifdef _OPENMP
        omp_init_lock(&mnode_lock);
#endif
    }

    Node(IndexType NewId )
        : BaseType()
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition()
    {
        KRATOS_ERROR <<  "Calling the default constructor for the node ... illegal operation!!" << std::endl;
        CreateSolutionStepData();

#ifdef _OPENMP
        omp_init_lock(&mnode_lock);
#endif
    }

    /// 1d constructor.
    Node(IndexType NewId, const DataType NewX)
        : BaseType(NewX)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition(NewX)
    {
        CreateSolutionStepData();

#ifdef _OPENMP
        omp_init_lock(&mnode_lock);
#endif
    }

    /// 2d constructor.
    Node(IndexType NewId, const DataType NewX, const DataType NewY)
        : BaseType(NewX, NewY)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition(NewX, NewY)
    {
        CreateSolutionStepData();

#ifdef _OPENMP
        omp_init_lock(&mnode_lock);
#endif
    }

    /// 3d constructor.
    Node(IndexType NewId, const DataType NewX, const DataType NewY, const DataType NewZ)
        : BaseType(NewX, NewY, NewZ)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition(NewX, NewY, NewZ)
    {
        CreateSolutionStepData();

#ifdef _OPENMP
        omp_init_lock(&mnode_lock);
#endif
    }

    /// Point constructor.
    Node(IndexType NewId, PointType const& rThisPoint)
        : BaseType(rThisPoint)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition(rThisPoint)
    {
        CreateSolutionStepData();

#ifdef _OPENMP
        omp_init_lock(&mnode_lock);
#endif
    }

    /** Copy constructor. Initialize this node with given node.*/
    Node(Node const& rOtherNode)
        : BaseType(rOtherNode)
        , IndexedObject(rOtherNode)
        , Flags(rOtherNode)
        , mData(rOtherNode.mData)
        , mSolutionStepsNodalData(rOtherNode.mSolutionStepsNodalData)
        , mInitialPosition(rOtherNode.mInitialPosition)
    {
        //TODO ... this copy constructor should be removed sometimes as it is often source of error
        //KRATOS_THROW_ERROR(std::logic_error, "copying Nodes is not allowed", "");

        // Deep copying the dofs
        for(typename DofsContainerType::const_iterator i_dof = rOtherNode.mDofs.begin() ; i_dof != rOtherNode.mDofs.end() ; i_dof++)
           pAddDof(*i_dof);

#ifdef _OPENMP
        omp_init_lock(&mnode_lock);
#endif
    }

    /** Copy constructor from a node with different dimension.*/
    template<SizeType TOtherDimension>
    Node(Node<TOtherDimension> const& rOtherNode)
        : BaseType(rOtherNode)
        , IndexedObject(rOtherNode)
        , Flags(rOtherNode)
        , mDofs(rOtherNode.mDofs)
        , mData(rOtherNode.mData)
        , mSolutionStepsNodalData(rOtherNode.mSolutionStepsNodalData)
        , mInitialPosition(rOtherNode.mInitialPosition)
    {

#ifdef _OPENMP
        omp_init_lock(&mnode_lock);
#endif
    }

    /** Copy constructor from a point with different dimension.*/
    template<SizeType TOtherDimension>
    Node(IndexType NewId, Point<TOtherDimension> const& rThisPoint)
        : BaseType(rThisPoint)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition(rThisPoint)
    {
        CreateSolutionStepData();

#ifdef _OPENMP
        omp_init_lock(&mnode_lock);
#endif
    }

    /**
     * Constructor using coordinates stored in given array. Initialize
    this point with the coordinates in the array. */
    template<class TVectorType>
    Node(IndexType NewId, vector_expression<TVectorType> const&  rOtherCoordinates)
        : BaseType(rOtherCoordinates)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition(rOtherCoordinates)
    {
        CreateSolutionStepData();

#ifdef _OPENMP
        omp_init_lock(&mnode_lock);
#endif
    }



    /** Constructor using coordinates stored in given std::vector. Initialize
    this point with the coordinates in the array. */
    Node(IndexType NewId, std::vector<DataType> const&  rOtherCoordinates)
        : BaseType(rOtherCoordinates)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition()
    {
        CreateSolutionStepData();

#ifdef _OPENMP
        omp_init_lock(&mnode_lock);
#endif
    }

    /// 3d with variables list and data constructor.
    Node(IndexType NewId, const DataType NewX, const DataType NewY, const DataType NewZ, VariablesList*  pVariablesList, BlockType const* ThisData, SizeType NewQueueSize = 1)
        : BaseType(NewX, NewY, NewZ)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData(pVariablesList,ThisData,NewQueueSize)
        , mInitialPosition(NewX, NewY, NewZ)
    {
#ifdef _OPENMP
        omp_init_lock(&mnode_lock);
#endif
    }


    /// Destructor.
    ~Node() override
    {
#ifdef _OPENMP
        omp_destroy_lock(&mnode_lock);
#endif
    }

    void SetId(IndexType NewId) override
    {
        IndexedObject::SetId(NewId);
        typename Node<TDimension>::DofsContainerType& my_dofs = (this)->GetDofs();
        for(typename Node<TDimension>::DofsContainerType::iterator iii = my_dofs.begin();    iii != my_dofs.end(); iii++)
        {
            iii->SetId(NewId);
        }
    }

#ifdef _OPENMP
    omp_lock_t& GetLock()
    {
        return mnode_lock;
    }
#endif

    inline void SetLock()
    {
        //does nothing if openMP is not present
#ifdef _OPENMP
        omp_set_lock(&mnode_lock);
#endif
    }

    inline void UnSetLock()
    {
        //does nothing if openMP is not present
#ifdef _OPENMP
        omp_unset_lock(&mnode_lock);
#endif
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Node& operator=(const Node& rOther)
    {
        BaseType::operator=(rOther);

    // Deep copying the dofs
        for(typename DofsContainerType::const_iterator i_dof = rOther.mDofs.begin() ; i_dof != rOther.mDofs.end() ; i_dof++)
           pAddDof(*i_dof);

        mData = rOther.mData;
        mSolutionStepsNodalData = rOther.mSolutionStepsNodalData;
        mInitialPosition = rOther.mInitialPosition;

        return *this;
    }

    /// Assignment operator.
    template<SizeType TOtherDimension>
    Node& operator=(const Node<TOtherDimension>& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator =(rOther);
        IndexedObject::operator=(rOther);
        for(typename DofsContainerType::const_iterator i_dof = rOther.mDofs.begin() ; i_dof != rOther.mDofs.end() ; i_dof++)
           pAddDof(*i_dof);

        mData = rOther.mData;
        mSolutionStepsNodalData = rOther.mSolutionStepsNodalData;
        mInitialPosition = rOther.mInitialPosition;

        return *this;
    }

    bool operator==(const Node& rOther)
    {
        return PointType::operator ==(rOther);
    }

    template<class TVariableType> typename TVariableType::Type& operator()(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        return GetSolutionStepValue(rThisVariable, SolutionStepIndex);
    }

    template<class TVariableType> typename TVariableType::Type& operator()(const TVariableType& rThisVariable)
    {
        return GetSolutionStepValue(rThisVariable);
    }

    template<class TDataType> TDataType& operator[](const Variable<TDataType>& rThisVariable)
    {
        return GetValue(rThisVariable);
    }

    template<class TDataType> const TDataType& operator[](const Variable<TDataType>& rThisVariable) const
    {
        return GetValue(rThisVariable);
    }

    template<class TAdaptorType> typename TAdaptorType::Type& operator[](const VariableComponent<TAdaptorType>& rThisVariable)
    {
        return GetValue(rThisVariable);
    }

    template<class TAdaptorType> const typename TAdaptorType::Type& operator[](const VariableComponent<TAdaptorType>& rThisVariable) const
    {
        return GetValue(rThisVariable);
    }

    DataType& operator[](IndexType ThisIndex)
    {
        return BaseType::operator[](ThisIndex);
    }

    DataType operator[](IndexType ThisIndex) const
    {
        return BaseType::operator[](ThisIndex);
    }

    ///@}
    ///@name Nodal Data
    ///@{

    void CreateSolutionStepData()
    {
        mSolutionStepsNodalData.PushFront();
//    VariablesListDataValueContainer temp(&Globals::DefaultVariablesList);
//    if(!mSolutionStepsNodalData.empty())
//        mSolutionStepsNodalData.push_front(temp);
//    else
//        mSolutionStepsNodalData.resize(1,temp);
//    mSolutionStepsNodalData->push_front(DataValueContainer());

    }

    void CloneSolutionStepData()
    {
        mSolutionStepsNodalData.CloneFront();
//    if(!mSolutionStepsNodalData.empty())
//        mSolutionStepsNodalData.push_front(mSolutionStepsNodalData[SourceSolutionStepIndex]);
    }

    void OverwriteSolutionStepData(IndexType SourceSolutionStepIndex, IndexType DestinationSourceSolutionStepIndex)
    {
        mSolutionStepsNodalData.AssignData(mSolutionStepsNodalData.Data(SourceSolutionStepIndex), DestinationSourceSolutionStepIndex);
    }

    void ClearSolutionStepsData()
    {
        mSolutionStepsNodalData.Clear();
    }

    void SetSolutionStepVariablesList(VariablesList* pVariablesList)
    {
        mSolutionStepsNodalData.SetVariablesList(pVariablesList);
    }

    VariablesListDataValueContainer& SolutionStepData()
    {
        return mSolutionStepsNodalData;
    }

    const VariablesListDataValueContainer& SolutionStepData() const
    {
        return mSolutionStepsNodalData;
    }

    DataValueContainer& Data()
    {
        return mData;
    }

    const DataValueContainer& Data() const
    {
        return mData;
    }

    DataValueContainer& GetData()
    {
        return mData;
    }

    const DataValueContainer& GetData() const
    {
        return mData;
    }

    template<class TVariableType> typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable)
    {
        return mSolutionStepsNodalData.GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable) const
    {
        return mSolutionStepsNodalData.GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        return mSolutionStepsNodalData.GetValue(rThisVariable, SolutionStepIndex);
    }

    template<class TVariableType> typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
    {
        return mSolutionStepsNodalData.GetValue(rThisVariable, SolutionStepIndex);
    }


    template<class TDataType> bool SolutionStepsDataHas(const Variable<TDataType>& rThisVariable) const
    {
        return mSolutionStepsNodalData.Has(rThisVariable);
    }
    template<class TAdaptorType> bool SolutionStepsDataHas(const VariableComponent<TAdaptorType>& rThisVariable) const
    {
        return mSolutionStepsNodalData.Has(rThisVariable);
    }

    //*******************************************************************************************
    //By Riccardo
    //very similar to the one before BUT throws an error if the variable does not exist
    template<class TVariableType> typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable)
    {
        return mSolutionStepsNodalData.FastGetValue(rThisVariable);
    }

    template<class TVariableType> const typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable) const
    {
        return mSolutionStepsNodalData.FastGetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        return mSolutionStepsNodalData.FastGetValue(rThisVariable, SolutionStepIndex);
    }

    template<class TVariableType> const typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
    {
        return mSolutionStepsNodalData.FastGetValue(rThisVariable, SolutionStepIndex);
    }

    template<class TVariableType> typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex, IndexType ThisPosition)
    {
        return mSolutionStepsNodalData.FastGetValue(rThisVariable, SolutionStepIndex, ThisPosition);
    }

    template<class TVariableType> typename TVariableType::Type& FastGetCurrentSolutionStepValue(const TVariableType& rThisVariable, IndexType ThisPosition)
    {
        return mSolutionStepsNodalData.FastGetCurrentValue(rThisVariable, ThisPosition);
    }
//*******************************************************************************************

    template<class TVariableType> typename TVariableType::Type& GetValue(const TVariableType& rThisVariable)
    {
        return mData.GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type const& GetValue(const TVariableType& rThisVariable) const
    {
        return mData.GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type& GetValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        if(!mData.Has(rThisVariable))
            return mSolutionStepsNodalData.GetValue(rThisVariable, SolutionStepIndex);

        return mData.GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type const& GetValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
    {
        if(!mData.Has(rThisVariable))
            return mSolutionStepsNodalData.GetValue(rThisVariable, SolutionStepIndex);

        return mData.GetValue(rThisVariable);
    }

    template<class TVariableType>
    void SetValue(const TVariableType& rThisVariable, typename TVariableType::Type const& rValue)
    {
        mData.SetValue(rThisVariable, rValue);
    }

    template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }
    template<class TAdaptorType> bool Has(const VariableComponent<TAdaptorType>& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }

    ///@}
    ///@name Operations
    ///@{

    template<class TVariableType>
    inline void Fix(const TVariableType& rDofVariable)
    {
        pAddDof(rDofVariable)->FixDof();
    }

    template<class TVariableType>
    inline void Free(const TVariableType& rDofVariable)
    {
        pAddDof(rDofVariable)->FreeDof();
    }

    IndexType GetBufferSize() const
    {
        return mSolutionStepsNodalData.QueueSize();
    }

    void SetBufferSize(IndexType NewBufferSize)
    {
        mSolutionStepsNodalData.Resize(NewBufferSize);
    }

    ///@}
    ///@name Access
    ///@{

    const PointType& GetInitialPosition() const
    {
        return mInitialPosition;
    }
    PointType& GetInitialPosition()
    {
        return mInitialPosition;
    }

    DataType& X0()
    {
        return mInitialPosition.X();
    }
    DataType& Y0()
    {
        return mInitialPosition.Y();
    }
    DataType& Z0()
    {
        return mInitialPosition.Z();
    }

    DataType X0() const
    {
        return mInitialPosition.X();
    }
    DataType Y0() const
    {
        return mInitialPosition.Y();
    }
    DataType Z0() const
    {
        return mInitialPosition.Z();
    }

    void SetInitialPosition(const PointType& NewInitialPosition)
    {
        mInitialPosition.X() = NewInitialPosition.X();
        mInitialPosition.Y() = NewInitialPosition.Y();
        mInitialPosition.Z() = NewInitialPosition.Z();
    }

    void SetInitialPosition(DataType X, DataType Y, DataType Z)
    {
        mInitialPosition.X() = X;
        mInitialPosition.Y() = Y;
        mInitialPosition.Z() = Z;
    }

    VariablesList* pGetVariablesList()
    {
        return mSolutionStepsNodalData.pGetVariablesList();
    }

    const VariablesList* pGetVariablesList() const
    {
        return mSolutionStepsNodalData.pGetVariablesList();
    }

    ///@}
    ///@name Dofs
    ///@{

    /**
     * @brief Get DoF  position with a given position
     * @param rDofVariable Name of the variable to search the position
     * @tparam TVariableType The variable type template argument
     * @return The position of the given DoF variable
     */
    template<class TVariableType>
    inline unsigned int GetDofPosition(TVariableType const& rDofVariable) const
    {
        typename DofsContainerType::iterator it=mDofs.find(rDofVariable.Key());
        return it - mDofs.begin();
    }

    /**
     * @brief Get dof with a given position. If not found it is search
     * @param rDofVariable Name of the variable
     * @param pos Position of the DoF
     * @tparam TVariableType The variable type template argument
     * @return The DoF associated to the given variable
     */
    template<class TVariableType>
    inline const DofType& GetDof(TVariableType const& rDofVariable, int pos) const
    {
        typename DofsContainerType::const_iterator it_begin = mDofs.begin();
        typename DofsContainerType::const_iterator it_end = mDofs.end();
        typename DofsContainerType::const_iterator it;
        //if the guess is exact return the guess
        if(pos < it_end-it_begin)
        {
            it = it_begin + pos;
            if( (it)->GetVariable() == rDofVariable)
                return *it;
        }

        //otherwise do a find
        it = mDofs.find(rDofVariable.Key());
        return *(it);
    }

    template<class TVariableType>
    inline DofType& GetDof(TVariableType const& rDofVariable, int pos)
    {
        typename DofsContainerType::iterator it_begin=mDofs.begin();
        typename DofsContainerType::iterator it_end=mDofs.end();
        typename DofsContainerType::iterator it;
        //if the guess is exact return the guess
        if(pos < it_end-it_begin)
        {
            it = it_begin + pos;
            if( (it)->GetVariable() == rDofVariable)
                return *it;
        }

        //otherwise do a find
        it = mDofs.find(rDofVariable.Key());
        return *(it);
    }

    /**
     * @brief Get DoF for a given variable
     * @param rDofVariable Name of the variable
     * @tparam TVariableType The variable type template argument
     * @return The DoF associated to the given variable
     */
    template<class TVariableType>
    inline const DofType& GetDof(TVariableType const& rDofVariable) const
    {
        typename DofsContainerType::const_iterator it=mDofs.find(rDofVariable.Key());
        if ( it!= mDofs.end() )
            return *it;

        KRATOS_ERROR << "Not existant DOF in node #" << Id() << " for variable : " << rDofVariable.Name();
    }

    /** retuns the Dof asociated with variable  */
    template<class TVariableType>
    inline DofType& GetDof(TVariableType const& rDofVariable)
    {
        typename DofsContainerType::iterator it=mDofs.find(rDofVariable.Key());
        if ( it!= mDofs.end() )
            return *it;

        KRATOS_ERROR << "Not existant DOF in node #" << Id() << " for variable : " << rDofVariable.Name();
    }

    /** retuns all of the Dofs  */
    DofsContainerType& GetDofs()
    {
        return mDofs;
    }

    /** retuns a counted pointer to the Dof asociated with variable  */
    template<class TVariableType>
    inline typename DofType::Pointer pGetDof(TVariableType const& rDofVariable)
    {
        typename DofsContainerType::iterator it=mDofs.find(rDofVariable.Key());
        if ( it!= mDofs.end() )
            return *(it.base());

        KRATOS_ERROR << "Not existant DOF in node #" << Id() << " for variable : " << rDofVariable.Name();
    }

    /**
     * @brief Get DoF counted pointer for a given variable
     * @param rDofVariable Name of the variable
     * @tparam TVariableType The variable type template argument
     * @return The DoF associated to the given variable
     */
    template<class TVariableType>
    inline const typename DofType::Pointer pGetDof(TVariableType const& rDofVariable) const
    {
        typename DofsContainerType::const_iterator it=mDofs.find(rDofVariable.Key());
        if ( it!= mDofs.end() )
            return *(it.base());

        KRATOS_ERROR <<  "Non-existent DOF in node #" << Id() << " for variable : " << rDofVariable.Name() << std::endl;
    }

    /**
     * @brief Get DoF counted pointer with a given position. If not found it is search
     * @param rDofVariable Name of the variable
     * @param Position Position of the DoF
     * @tparam TVariableType The variable type template argument
     * @return The DoF associated to the given variable
     */
    template<class TVariableType>
    inline typename DofType::Pointer pGetDof(
        TVariableType const& rDofVariable,
        int Position
        ) const
    {
        const auto it_begin = mDofs.begin();
        const auto it_end = mDofs.end();
        // If the guess is exact return the guess
        if(Position < it_end-it_begin) {
            auto it_dof = it_begin + Position;
            if( (*it_dof)->GetVariable() == rDofVariable) {
                return (*it_dof).get();
            }
        }

        // Otherwise do a find
        for(auto it_dof = it_begin; it_dof != it_end; ++it_dof){
            if((*it_dof)->GetVariable() == rDofVariable){
                return (*it_dof).get();
            }
        }

        KRATOS_ERROR <<  "Non-existent DOF in node #" << Id() << " for variable : " << rDofVariable.Name() << std::endl;
    }

    /** adds a Dof to the node and return new added dof or existed one. */
    template<class TVariableType>
    inline typename DofType::Pointer pAddDof(TVariableType const& rDofVariable)
    {
        KRATOS_TRY_LEVEL_3

        typename DofsContainerType::iterator i_dof = mDofs.find(rDofVariable);
        if(i_dof != mDofs.end())
            return *(i_dof.base());
        return *(mDofs.insert(mDofs.begin(), DofType(Id(), &mSolutionStepsNodalData, rDofVariable)).base());

        KRATOS_CATCH_LEVEL_3(*this);

    }

    /** adds a Dof to the node and return new added dof or existed one. */
    inline typename DofType::Pointer pAddDof(DofType const& SourceDof)
    {
        KRATOS_TRY_LEVEL_3

        typename DofsContainerType::iterator i_dof = mDofs.find(SourceDof.GetVariable());
        if(i_dof != mDofs.end())
            return *(i_dof.base());

        typename DofType::Pointer p_new_dof(  *(mDofs.insert(mDofs.begin(), SourceDof).base())  );

        p_new_dof->SetId(Id());
//      p_new_dof->Id() = Id();

        p_new_dof->SetSolutionStepsData(&mSolutionStepsNodalData);

        return p_new_dof;

        KRATOS_CATCH_LEVEL_3(*this);
    }

    /** adds a Dof to the node and return new added dof or existed one. */
    template<class TVariableType, class TReactionType>
    inline typename DofType::Pointer pAddDof(TVariableType const& rDofVariable, TReactionType const& rDofReaction)
    {
        KRATOS_TRY_LEVEL_3

        typename DofsContainerType::iterator i_dof = mDofs.find(rDofVariable);
        if(i_dof != mDofs.end())
        {
            i_dof->SetReaction(rDofReaction);
            return *(i_dof.base());
        }
        return *(mDofs.insert(mDofs.begin(), DofType(Id(), &mSolutionStepsNodalData, rDofVariable, rDofReaction)).base());

        KRATOS_CATCH_LEVEL_3(*this);
    }

    /** adds a Dof to the node and return new added dof or existed one. */
    template<class TVariableType>
    inline DofType& AddDof(TVariableType const& rDofVariable)
    {
        KRATOS_TRY_LEVEL_3

        typename DofsContainerType::iterator i_dof = mDofs.find(rDofVariable);
        if(i_dof != mDofs.end())
            return *i_dof;
        return *(mDofs.insert(mDofs.begin(), DofType(Id(), &mSolutionStepsNodalData, rDofVariable)));

        KRATOS_CATCH_LEVEL_3(*this);
    }

    /** adds a Dof to the node and return new added dof or existed one. */
    template<class TVariableType, class TReactionType>
    inline DofType& AddDof(TVariableType const& rDofVariable, TReactionType const& rDofReaction)
    {
        KRATOS_TRY_LEVEL_3

        typename DofsContainerType::iterator i_dof = mDofs.find(rDofVariable);
        if(i_dof != mDofs.end())
            return *i_dof;
        return *(mDofs.insert(mDofs.begin(), DofType(Id(), &mSolutionStepsNodalData, rDofVariable, rDofReaction)));

        KRATOS_CATCH_LEVEL_3(*this);
    }

    ///@}
    ///@name Inquiry
    ///@{

    /** Return true if the dof of freedom is present on the node */
    inline bool HasDofFor(const VariableData& rDofVariable) const
    {
        return (mDofs.find(rDofVariable) != mDofs.end());
    }

    inline  bool IsFixed(const VariableData& rDofVariable) const
    {
        typename DofsContainerType::const_iterator i;
        return (((i= mDofs.find(rDofVariable)) == mDofs.end()) ? false : i->IsFixed());
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Node #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
        if(!mDofs.empty())
            rOStream << std::endl << "    Dofs :" << std::endl;
        for(typename DofsContainerType::const_iterator i = mDofs.begin() ; i != mDofs.end() ; i++)
            rOStream << "        " << i->Info() << std::endl;
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

#ifdef _OPENMP
    omp_lock_t mnode_lock;
#endif


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

    /** storage for the dof of the node */
    DofsContainerType  mDofs;

    /** A pointer to data related to this node. */
    DataValueContainer mData;

    SolutionStepsNodalDataContainerType mSolutionStepsNodalData;

    //SolutionStepsNeighboursContainerType mSolutionStepsNeighbours;

    ///Initial Position of the node
    PointType mInitialPosition;


    ///@}
    ///@name Private Operators
    ///@{



    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
//    int size = rSerializer.GetBuffer().end() - rSerializer.GetBuffer().begin();
//    KRATOS_WATCH(rSerializer.GetBuffer().end() - rSerializer.GetBuffer().begin());
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PointType );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
        rSerializer.save("Data", mData);
        const SolutionStepsNodalDataContainerType* p_solution_steps_nodal_data = &mSolutionStepsNodalData;
        // I'm saving it as pointer so the dofs pointers will point to it as stored pointer. Pooyan.
        rSerializer.save("Solution Steps Nodal Data", p_solution_steps_nodal_data);
        rSerializer.save("Initial Position", mInitialPosition);
        rSerializer.save("Data", mDofs);

//    KRATOS_WATCH((rSerializer.GetBuffer().end() - rSerializer.GetBuffer().begin())-size);
    }

    void load(Serializer& rSerializer) override
    {
//    int size = rSerializer.GetBuffer().end() - rSerializer.GetBuffer().begin();
//    KRATOS_WATCH(rSerializer.GetBuffer().end() - rSerializer.GetBuffer().begin());
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PointType );
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject );
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
        rSerializer.load("Data", mData);
        SolutionStepsNodalDataContainerType* p_solution_steps_nodal_data = &mSolutionStepsNodalData;
        rSerializer.load("Solution Steps Nodal Data", p_solution_steps_nodal_data);
        rSerializer.load("Initial Position", mInitialPosition);
        rSerializer.load("Data", mDofs);
//    KRATOS_WATCH(size- (rSerializer.GetBuffer().end() - rSerializer.GetBuffer().begin()));
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

}; // Class Node

///@}

// template class KRATOS_API(KRATOS_CORE) KratosComponents<Node<3, KRATOS_DOUBLE_TYPE> >;
// template class KRATOS_API(KRATOS_CORE) KratosComponents<Node<3, KRATOS_DOUBLE_TYPE>::Pointer >;

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<std::size_t TDimension, class TDofType>
inline std::istream& operator >> (std::istream& rIStream,
                                  Node<TDimension, TDofType>& rThis);

/// output stream function
template<std::size_t TDimension, class TDofType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Node<TDimension, TDofType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : ";
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

//*********************************************************************************
//*********************************************************************************
//*********************************************************************************
//definition of the NEIGHBOUR_NODES variable
//*********************************************************************************
//*********************************************************************************
//*********************************************************************************

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

KRATOS_DEFINE_VARIABLE(WeakPointerVector<Node<3> >, NEIGHBOUR_NODES)
KRATOS_DEFINE_VARIABLE(WeakPointerVector<Node<3> >, FATHER_NODES)

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT


//     namespace Globals
//     {
//  extern Node<3> DefaultNode3;
//     }

}  // namespace Kratos.

#endif // KRATOS_NODE_H_INCLUDED  defined
