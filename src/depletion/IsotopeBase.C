#include "IsotopeBase.h"

#include <cmath>
#include <set>

IsotopeBase::IsotopeBase(unsigned id, pugi::xml_node isotope_entry)
 : _id(id)
{
    this->unpackXMLNode(isotope_entry);

    //Including self, dont care about double counting rows
    //even if its less efficient
    PetscInt row_ids[this->_children.size() + 1];
    this->_row_ids = row_ids;

    Real mat_entries[this->_children.size() + 1];
    this->_mat_entries = mat_entries;

    this->_row_size = 1; //self
}

void IsotopeBase::unpackXMLNode(pugi::xml_node node)
{
//reading xml node

    this->_name = node.attribute("name").as_string();


    //TRANSMUTE
    //PetscInt num_nx = node.attribute("reactions").as_int();

    // not yet sure how to do xs with njoy, might switch this to std::vector<enum>
    std::vector<std::string> nx_types;

    for (auto child : node.children("reaction"))
    {
        this->_nx_children.push_back( child.attribute("target").as_string() );
        this->_nx_types.push_back( child.attribute( "type" ).as_string() );
    }
    

    //DECAY
    //PetscInt num_dec = node.attribute("decay_modes").as_int();

    Real decay_constant = std::log(2) / node.attribute("half_life").as_double();

    for (auto child : node.children("decay"))
    {
        this->_decay_children.push_back( child.attribute( "target" ).as_string() );
        this->_decay_entries.push_back( decay_constant * 
                child.attribute( "branching_ratio" ).as_double() );
    }
 
//unpacking time

    this->_children = this->_decay_children;
    this->_children.insert(this->_children.end(), 
            this->_nx_children.begin(), this->_nx_children.end());
}


void IsotopeBase::setChildrenIDs(PetscInt ids[])
{
    this->_row_ids[0] = this->_id;

    for (int i = 0; i < this->_children.size(); i++)
    {
        this->_row_ids[i + 1] = ids[i];
    }
}

void
IsotopeBase::setMatrixEntries()
{
    
}