#pragma once

#include <vector>
#include <tuple>
#include <string>

#include "libmesh/libmesh_common.h"
#include "pugixml.h"

#include "MooseVariableBase.h"
#include "MooseError.h"


class IsotopeBase
{
public:
    IsotopeBase(unsigned id,
               pugi::xml_node isotope_entry);

    void unpackXMLNode(pugi::xml_node isotope_entry);

    std::vector<std::string> getChildrenNames() const {return this->_children;}

    void setChildrenIDs(PetscInt ids[]);

/**
 * GETTER FUNCTIONS
 *  - Called by DepletionMatrix::evaluateMatrix() 
 *    or DepletionMatrix::readChainXML()
 */   

    /**
     * @return Isotope's id
     * @paragraph Everytime this function is called, a parent is accesing
     * this isotope's ID, and so the row size of this Isotope increases by 1
     */
    PetscInt getID() {this->_row_size += 1; return this->_id;}

    /**
     * @return Isotope's name. "Xe135"
     */
    std::string getName() const {return this->_name;}

    /**
     * @return number of rows this Isotope touches
     * 
     */
    PetscInt getNumRows() const {return this->_children.size() + 1;}

    /**
     * @return number of columns this Isotope touches
     */
    PetscInt getNumCols() const {return 1;}

    /**
     * @return an int[] of the index of each row this isotope touches
     * Each entry in the returned array relates to the corresponding row
     * index returned by getRowIDs()
     */
    PetscInt* getRowIDs() const {return this->_row_ids;}

    /**
     * @return an int[] of the index of each column this isotope touches
     * Each entry in the returned array relates to the corresponding row
     * index returned by getRowIDs()
     */
    PetscInt* getColIDs() {return &this->_id;}

    /**
     * @return an Real[] of the values for each (i,j) in (row_ids, col_ids)
     * entries returned include:
     *  (i)   the loss term from this isotope due to decay
     *  (ii)  the source term for children from this isotope's decay
     *  (iii) the source term for children from this isotope's (n,x) reactions
     */
    Real* getMatrixEntries() const {return this->_mat_entries;}

    /**
     * @return the number of entries in this isotope's row. 
     * This number is the number of parents + 1
     */
    PetscInt getRowSize() const {return this->_row_size;}
    
protected:

    virtual void setMatrixEntries();

    PetscInt _id;

    std::string _name;

    std::vector<std::string> _nx_children;
    std::vector<std::string> _nx_types;
    std::vector<Real> _nx_entries;

    std::vector<std::string> _decay_children;
    std::vector<Real> _decay_entries;

    std::vector<std::string> _children;

    PetscInt* _row_ids;
    Real* _mat_entries;

    PetscInt _row_size;
};
