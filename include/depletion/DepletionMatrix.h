#pragma once

#include "ThreadedGeneralUserObject.h"
#include "MooseVariableBase.h"
#include "IsotopeBase.h"
#include "json.h"
#include "pugixml.h"

#include <petscmat.h>

class DepletionMatrix : public ThreadedGeneralUserObject
{
public:

    static InputParameters validParams();

    DepletionMatrix(const InputParameters & parameters);

    virtual void execute() override;
    
    virtual void initialize() final {};

    virtual void finalize() final {};

    Mat getDepletionMatrix() const {return this->_depl_mat;}

    std::vector<std::string> getIsotopeNames()
        const {return this->_isotope_names;}

protected:

    /**
     * Create the NNX for the depletion matrix.
     * This tells PetSC how much memory to allocate
     * for the matrix.
     */
    void setNNZ();

    /**
     * Sets up all Isotopes with their data from the 
     * given decay chain file
     * @param file the depletion chain file to read
     * 
     */
    void readChainXML(std::string file);

    /**
     * Fills in the Depletion matrix with data from 
     * the Isotopes.
     */
    void evaluateMatrix();

    std::vector<IsotopeBase> _isotopes;

    std::vector<std::string> _isotope_names;

    PetscInt _num_isos;

    PetscInt * _nnz;

    Mat _depl_mat;
};
