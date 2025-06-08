#include "DepletionMatrix.h"

registerMooseObject("DittoApp", DepletionMatrix);

InputParameters
DepletionMatrix::validParams()
{
    auto params = ThreadedGeneralUserObject::validParams();
    params.addClassDescription("This class implements the depletion matrix"
        "used by the depletion solvers.");
    params.addRequiredParam<std::string>("isotope_file", 
        "The file holding all of the decay and transmutation data for the simulation");
    params.addRequiredParam<bool>("evaluate_xs",
        "Whether to evaluate Isotope cross sections on the fly using NJOY.");
    return params;
}

DepletionMatrix::DepletionMatrix(const InputParameters & parameters)
: ThreadedGeneralUserObject(parameters)
{
    bool eval_xs = getParam<bool>("evaluate_xs");
    std::string file = getParam<std::string>("isotope_file");
    
    readChainXML(file, eval_xs)

    _num_isos = _isotopes.size();
    
    setNNZ();
    // Creates sparse PETSc matrix (SeqAIJ). Allocates the memory of the matrix with nnz not nz
    PetscCall(MatCreateSeqAIJ(PETSC_COMM_SELF, _num_isos, _num_isos, 0, this->_nnz, &_depl_matrix));

    evaluateMatrix();
}

void DepletionMatrix::execute()
{
    evaluateMatrix();
}

void
DepletionMatrix::setNNZ()
{
    PetscInt nnz[_num_isos];

    for (int i = 0; i < _num_isos; i++)
    {
        nnz[i] = _isotopes.at(i).getRowSize();
    }

    this->_nnz = nnz;
}

void
DepletionMatrix::evaluateMatrix()
{
    PetscCall( MatAssemblyBegin(_depl_matrix, MAT_FINAL_ASSEMBLY) );

    PetscCall( MatZeroEntries(_depl_matrix) );
    
    for (auto isotope : _isotopes)
    {
        PetscCall( MatSetValues(/*Mat*/ this->_depl_matrix, 
                                /*Num Rows*/ isotope.getNumRows(),
                                /*Row Indices*/ isotope.getRowIDs(), //this should be arr of len num cols
                                /*Num Cols*/ isotope.getNumCols(),
                                /*Col Indices*/ isotope.getColIDs(),
                                /*Values*/ isotope.getMatrixEntries(),
                                /*Insert Mode*/ INSERT_VALUES) );
    }

    PetscCall( MatAssemblyEnd(_depl_matrix, MAT_FINAL_ASSEMBLY) );
}

void 
DepletionMatrix::readChainXML(std::string file, bool eval_xs)
{
    pugi::xml_document doc;
    if (!doc.load_file(file))
    {
    mooseError("Failed to open "+file+", ensure it is formatted properly!");
    }

    unsigned id = 0;
    for (auto & node : doc.child("depletion_chain"))
    {
        auto iso = IsotopeBase(id, makeNewVariable(), node);
        this->_isotopes.push_back( iso );
        this->_isotope_names.push_back( iso.getName() );
        id += 1;
    }
}
