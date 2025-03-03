//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_FIND_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED )
#define  KRATOS_FIND_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"


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
class FindElementalNeighboursProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindElementalNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindElementalNeighboursProcess);

    typedef ModelPart::NodesContainerType NodesContainerType;
    typedef ModelPart::ElementsContainerType ElementsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// avg_elems ------ expected number of neighbour elements per node.,
    /// avg_nodes ------ expected number of neighbour Nodes
    /// the better the guess for the quantities above the less memory occupied and the fastest the algorithm
    FindElementalNeighboursProcess(ModelPart& model_part, int TDim, unsigned int avg_elems = 10)
        : mr_model_part(model_part)
    {
        mavg_elems = avg_elems;
        mTDim=TDim;
// 	mavg_nodes = avg_nodes;
    }

    /// Destructor.
    ~FindElementalNeighboursProcess() override
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        NodesContainerType& rNodes = mr_model_part.Nodes();
        ElementsContainerType& rElems = mr_model_part.Elements();

        //first of all the neighbour nodes and elements array are initialized to the guessed size
        //and empties the old entries
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            (in->GetValue(NEIGHBOUR_ELEMENTS)).reserve(mavg_elems);
            WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(),rE.end() );
        }
        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
        {
            if(!ie->Is(ACTIVE)) continue;
            (ie->GetValue(NEIGHBOUR_ELEMENTS)).reserve(NumberOfNeighbours(ie->GetGeometry().GetGeometryType()));
            WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(),rE.end() );
        }

        //add the neighbour elements to all the nodes in the mesh
        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
        {
            if(!ie->Is(ACTIVE)) continue;
            Element::GeometryType& pGeom = ie->GetGeometry();
            for(unsigned int i = 0; i < pGeom.size(); i++)
            {
                //KRATOS_WATCH( pGeom[i] );
                (pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( *(ie.base()) ) );
                //KRATOS_WATCH( (pGeom[i].GetValue(NEIGHBOUR_CONDITIONS)).size() );
            }
        }

        //adding the neighbouring conditions to the condition
        //loop over faces
        if (mTDim==2)
        {
            for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
            {
                if(!ie->Is(ACTIVE)) continue;

                //face nodes
                Geometry<Node<3> >& geom = (ie)->GetGeometry();

                if(geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3 || geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D6)
                {
                    //vector of the 3 faces around the given face
                    (ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(3);
                    WeakPointerVector< Element >& neighb_elems = ie->GetValue(NEIGHBOUR_ELEMENTS);
                    //neighb_face is the vector containing pointers to the three faces around ic
                    //neighb_face[0] = neighbour face over edge 1-2 of element ic;
                    //neighb_face[1] = neighbour face over edge 2-0 of element ic;
                    //neighb_face[2] = neighbour face over edge 0-1 of element ic;
                    neighb_elems(0) = CheckForNeighbourElems(geom[1].Id(), geom[2].Id(), geom[1].GetValue(NEIGHBOUR_ELEMENTS), ie);
                    neighb_elems(1) = CheckForNeighbourElems(geom[2].Id(), geom[0].Id(), geom[2].GetValue(NEIGHBOUR_ELEMENTS), ie);
                    neighb_elems(2) = CheckForNeighbourElems(geom[0].Id(), geom[1].Id(), geom[0].GetValue(NEIGHBOUR_ELEMENTS), ie);
                }
                else if(geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4
                    || geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8
                    || geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9)
                {
                    //vector of the 3 faces around the given face
                    (ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(4);
                    WeakPointerVector< Element >& neighb_elems = ie->GetValue(NEIGHBOUR_ELEMENTS);
                    //neighb_face is the vector containing pointers to the three faces around ic
                    //neighb_face[0] = neighbour face over edge 1-2 of element ic;
                    //neighb_face[1] = neighbour face over edge 2-3 of element ic;
                    //neighb_face[2] = neighbour face over edge 3-0 of element ic;
                    //neighb_face[3] = neighbour face over edge 0-1 of element ic;
                    neighb_elems(0) = CheckForNeighbourElems(geom[1].Id(), geom[2].Id(), geom[1].GetValue(NEIGHBOUR_ELEMENTS), ie);
                    neighb_elems(1) = CheckForNeighbourElems(geom[2].Id(), geom[3].Id(), geom[2].GetValue(NEIGHBOUR_ELEMENTS), ie);
                    neighb_elems(2) = CheckForNeighbourElems(geom[3].Id(), geom[0].Id(), geom[3].GetValue(NEIGHBOUR_ELEMENTS), ie);
                    neighb_elems(3) = CheckForNeighbourElems(geom[0].Id(), geom[1].Id(), geom[0].GetValue(NEIGHBOUR_ELEMENTS), ie);
                }
                else
                {
                    std::stringstream ss;
                    ss << "The geometry type " << static_cast<int>(geom.GetGeometryType()) << " is invalid";
                    KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
                }
            }
        }
        if (mTDim==3)
        {
            for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
            {
                if(!ie->Is(ACTIVE)) continue;

                //face nodes
                Geometry<Node<3> >& geom = (ie)->GetGeometry();

                if(geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4
                    || geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10)
                {
                    //vector of the 3 faces around the given face
                    (ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(4);
                    WeakPointerVector< Element >& neighb_elems = ie->GetValue(NEIGHBOUR_ELEMENTS);
                    //neighb_face is the vector containing pointers to the three faces around ic
                    //neighb_face[0] = neighbour face over edge 1-2 of element ic;
                    //neighb_face[1] = neighbour face over edge 2-0 of element ic;
                    //neighb_face[2] = neighbour face over edge 0-1 of element ic;
                    neighb_elems(0) = CheckForNeighbourElemsTetra(geom[1].Id(), geom[2].Id(), geom[3].Id(), geom[1].GetValue(NEIGHBOUR_ELEMENTS), ie);
                    neighb_elems(1) = CheckForNeighbourElemsTetra(geom[2].Id(), geom[3].Id(), geom[0].Id(), geom[2].GetValue(NEIGHBOUR_ELEMENTS), ie);
                    neighb_elems(2) = CheckForNeighbourElemsTetra(geom[3].Id(), geom[0].Id(), geom[1].Id(), geom[3].GetValue(NEIGHBOUR_ELEMENTS), ie);
                    neighb_elems(3) = CheckForNeighbourElemsTetra(geom[0].Id(), geom[1].Id(), geom[2].Id(), geom[0].GetValue(NEIGHBOUR_ELEMENTS), ie);
                }
                else if(geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8
                    || geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D20
                    || geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D27)
                {
                    //vector of the 3 faces around the given face
                    (ie->GetValue(NEIGHBOUR_ELEMENTS)).resize(6);
                    WeakPointerVector< Element >& neighb_elems = ie->GetValue(NEIGHBOUR_ELEMENTS);
                    //neighb_face is the vector containing pointers to the three faces around ic
                    neighb_elems(0) = CheckForNeighbourElemsHexa(geom[3].Id(), geom[2].Id(), geom[1].Id(), geom[0].Id(), geom[3].GetValue(NEIGHBOUR_ELEMENTS), ie);
                    neighb_elems(1) = CheckForNeighbourElemsHexa(geom[0].Id(), geom[1].Id(), geom[5].Id(), geom[4].Id(), geom[0].GetValue(NEIGHBOUR_ELEMENTS), ie);
                    neighb_elems(2) = CheckForNeighbourElemsHexa(geom[2].Id(), geom[6].Id(), geom[5].Id(), geom[1].Id(), geom[2].GetValue(NEIGHBOUR_ELEMENTS), ie);
                    neighb_elems(3) = CheckForNeighbourElemsHexa(geom[7].Id(), geom[6].Id(), geom[2].Id(), geom[3].Id(), geom[7].GetValue(NEIGHBOUR_ELEMENTS), ie);
                    neighb_elems(4) = CheckForNeighbourElemsHexa(geom[7].Id(), geom[3].Id(), geom[0].Id(), geom[4].Id(), geom[7].GetValue(NEIGHBOUR_ELEMENTS), ie);
                    neighb_elems(5) = CheckForNeighbourElemsHexa(geom[4].Id(), geom[5].Id(), geom[6].Id(), geom[7].Id(), geom[4].GetValue(NEIGHBOUR_ELEMENTS), ie);
                }
                else
                {
                    std::stringstream ss;
                    ss << "The geometry type " << static_cast<int>(geom.GetGeometryType()) << " is invalid";
                    KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
                }
            }
        }
    }


    void ClearNeighbours()
    {
        NodesContainerType& rNodes = mr_model_part.Nodes();
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(),rE.end());
        }
        ElementsContainerType& rElems = mr_model_part.Elements();
        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
        {
            WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(),rE.end());
        }

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
        return "FindElementalNeighboursProcess";
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
    ModelPart& mr_model_part;
    unsigned int mavg_elems;
    int mTDim;
// 	unsigned int mavg_nodes;


    ///@}
    ///@name Private Operators
    ///@{

    //******************************************************************************************
    //******************************************************************************************
    template< class TDataType > void  AddUniqueWeakPointer
    (WeakPointerVector< TDataType >& v, const typename TDataType::WeakPointer candidate)
    {
        typename WeakPointerVector< TDataType >::iterator i = v.begin();
        typename WeakPointerVector< TDataType >::iterator endit = v.end();
        while ( i != endit && (i)->Id() != (candidate.lock())->Id())
        {
            i++;
        }
        if( i == endit )
        {
            v.push_back(candidate);
        }

    }

    Element::WeakPointer CheckForNeighbourElems (const unsigned int& Id_1, const unsigned int& Id_2,
        WeakPointerVector< Element >& neighbour_elem, ElementsContainerType::iterator elem)
    {
        //look for the faces around node Id_1
        for( WeakPointerVector< Element >::iterator i = neighbour_elem.begin(); i != neighbour_elem.end(); i++)
        {
            //look for the nodes of the neighbour faces
            Geometry<Node<3> >& neigh_elem_geometry = (i)->GetGeometry();
            for( unsigned int node_i = 0 ; node_i < neigh_elem_geometry.size(); node_i++)
            {
                if (neigh_elem_geometry[node_i].Id() == Id_2)
                {
                    if(i->Id() != elem->Id())
                    {
                        return *(i.base());
                    }
                }
            }
        }
        return *(elem.base());
    }

    Element::WeakPointer CheckForNeighbourElemsTetra (const unsigned int& Id_1,
        const unsigned int& Id_2, const unsigned int& Id_3,
        WeakPointerVector< Element >& neighbour_elem, ElementsContainerType::iterator elem)
    {
        //look for the faces around node Id_1
        for( WeakPointerVector< Element >::iterator i = neighbour_elem.begin(); i != neighbour_elem.end(); i++)
        {
            //look for the nodes of the neighbour faces
            Geometry<Node<3> >& neigh_elem_geometry = (i)->GetGeometry();
            for( unsigned int node_i = 0 ; node_i < neigh_elem_geometry.size(); node_i++)
            {
                if (neigh_elem_geometry[node_i].Id() == Id_2)
                {
                    for( unsigned int node_j = 0 ; node_j < neigh_elem_geometry.size(); node_j++)
                    {
                        if (neigh_elem_geometry[node_j].Id() == Id_3)
                            if(i->Id() != elem->Id())
                            {
                                return *(i.base());
                            }
                    }
                }
            }
        }
        return *(elem.base());
    }

    Element::WeakPointer CheckForNeighbourElemsHexa (const unsigned int& Id_1,
        const unsigned int& Id_2, const unsigned int& Id_3, const unsigned int& Id_4,
        WeakPointerVector< Element >& neighbour_elem, ElementsContainerType::iterator elem)
    {
        //look for the faces around node Id_1
        for( WeakPointerVector< Element >::iterator i = neighbour_elem.begin(); i != neighbour_elem.end(); i++)
        {
            //look for the nodes of the neighbour faces
            Geometry<Node<3> >& neigh_elem_geometry = (i)->GetGeometry();
            for( unsigned int node_i = 0 ; node_i < neigh_elem_geometry.size(); node_i++)
            {
                if (neigh_elem_geometry[node_i].Id() == Id_2)
                {
                    for( unsigned int node_j = 0 ; node_j < neigh_elem_geometry.size(); node_j++)
                    {
                        if (neigh_elem_geometry[node_j].Id() == Id_3)
                        {
                            for( unsigned int node_k = 0 ; node_k < neigh_elem_geometry.size(); node_k++)
                            {
                                if (neigh_elem_geometry[node_k].Id() == Id_4)
                                    if(i->Id() != elem->Id())
                                    {
                                        return *(i.base());
                                    }
                            }
                        }
                    }
                }
            }
        }
        return *(elem.base());
    }

    /// Get the number of neighbours of a specific geometry type. It is typically the same as number of element edges.
    static int NumberOfNeighbours(const GeometryData::KratosGeometryType& rGeometryType)
    {
        switch(rGeometryType)
        {
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3: return 3;
            case GeometryData::KratosGeometryType::Kratos_Triangle2D6: return 3;
            case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4: return 4;
            case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8: return 4;
            case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9: return 4;
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4: return 4;
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10: return 4;
            case GeometryData::KratosGeometryType::Kratos_Hexahedra3D8: return 6;
            case GeometryData::KratosGeometryType::Kratos_Hexahedra3D20: return 6;
            case GeometryData::KratosGeometryType::Kratos_Hexahedra3D27: return 6;
            default: return 0;
        }
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
    FindElementalNeighboursProcess& operator=(FindElementalNeighboursProcess const& rOther);

    /// Copy constructor.
    //FindElementalNeighboursProcess(FindConditionsNeighboursProcess const& rOther);


    ///@}

}; // Class FindConditionsNeighboursProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_FIND_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED  defined
