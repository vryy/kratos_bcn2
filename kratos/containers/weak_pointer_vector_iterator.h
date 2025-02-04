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

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_WEAK_POINTER_VECTOR_ITERATOR_H_INCLUDED )
#define  KRATOS_WEAK_POINTER_VECTOR_ITERATOR_H_INCLUDED



// System includes


// External includes
#include <boost/iterator/iterator_adaptor.hpp>


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

/// Short class definition.
/** Detail class definition.
*/
template<class TIteratorType, class TDataType>
class WeakPointerVectorIterator
    : public boost::iterator_adaptor<WeakPointerVectorIterator<TIteratorType, TDataType>,
      TIteratorType, TDataType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of WeakPointerVectorIterator
    KRATOS_CLASS_POINTER_DEFINITION(WeakPointerVectorIterator);

    typedef boost::iterator_adaptor<WeakPointerVectorIterator,
            TIteratorType, TDataType> BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    WeakPointerVectorIterator() {}

    WeakPointerVectorIterator(TIteratorType NewIterator) :BaseType(NewIterator) {}

    WeakPointerVectorIterator(WeakPointerVectorIterator const & NewIterator) :BaseType(NewIterator.base()) {}

    //template<class TOtherIteratorType>
    //WeakPointerVectorIterator(WeakPointerVectorIterator<TIteratorType> const & NewIterator) :BaseType(NewIterator.base()) {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    friend class boost::iterator_core_access;
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

    typename BaseType::reference dereference() const
    {
        return *(boost::shared_ptr<TDataType>(*(this->base())));
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

}; // Class WeakPointerVectorIterator

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_WEAK_POINTER_VECTOR_ITERATOR_H_INCLUDED  defined 


