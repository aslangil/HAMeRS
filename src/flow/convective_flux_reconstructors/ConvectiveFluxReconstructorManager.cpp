#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructorManager.hpp"

ConvectiveFluxReconstructorManager::ConvectiveFluxReconstructorManager(
    const std::string& object_name,
    const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db,
    const std::string& convective_flux_reconstructor_str):
        d_object_name(object_name),
        d_convective_flux_reconstructor_db(convective_flux_reconstructor_db)
{
    TBOX_ASSERT(!object_name.empty());
    
    if (convective_flux_reconstructor_str == "FIRST_ORDER_LLF")
    {
        d_convective_flux_reconstructor_type = CONVECTIVE_FLUX_RECONSTRUCTOR::FIRST_ORDER_LLF;
    }
    else if (convective_flux_reconstructor_str == "FIRST_ORDER_HLLC")
    {
        d_convective_flux_reconstructor_type = CONVECTIVE_FLUX_RECONSTRUCTOR::FIRST_ORDER_HLLC;
    }
    else if (convective_flux_reconstructor_str == "WCNS5_JS_HLLC_HLL")
    {
        d_convective_flux_reconstructor_type = CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS5_JS_HLLC_HLL;
    }
    else if (convective_flux_reconstructor_str == "WCNS5_Z_HLLC_HLL")
    {
        d_convective_flux_reconstructor_type = CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS5_Z_HLLC_HLL;
    }
    else if (convective_flux_reconstructor_str == "WCNS6_CU_M2_HLLC_HLL")
    {
        d_convective_flux_reconstructor_type = CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS6_CU_M2_HLLC_HLL;
    }
    else if (convective_flux_reconstructor_str == "WCNS6_LD_HLLC_HLL")
    {
        d_convective_flux_reconstructor_type = CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS6_LD_HLLC_HLL;
    }
    else if (convective_flux_reconstructor_str == "WCNS6_TEST")
    {
        d_convective_flux_reconstructor_type = CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS6_TEST;
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Unknown convective_flux_reconstructor string = '"
            << convective_flux_reconstructor_str
            << "' found in input."
            << std::endl);        
    }
}


/*
 * Create the convective flux reconstructor.
 */
boost::shared_ptr<ConvectiveFluxReconstructor>
ConvectiveFluxReconstructorManager::createConvectiveFluxReconstructor(
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const boost::shared_ptr<FlowModel>& flow_model)
{
    boost::shared_ptr<ConvectiveFluxReconstructor> conv_flux_reconstructor;
    
    switch (d_convective_flux_reconstructor_type)
    {
        case CONVECTIVE_FLUX_RECONSTRUCTOR::FIRST_ORDER_LLF:
        {
            conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorFirstOrderLLF(
                "d_convective_flux_reconstructor",
                dim,
                grid_geometry,
                flow_model->getNumberOfEquations(),
                num_species,
                flow_model,
                d_convective_flux_reconstructor_db));
            
            break;
        }
        case CONVECTIVE_FLUX_RECONSTRUCTOR::FIRST_ORDER_HLLC:
        {
            conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorFirstOrderHLLC(
                "d_convective_flux_reconstructor",
                dim,
                grid_geometry,
                flow_model->getNumberOfEquations(),
                num_species,
                flow_model,
                d_convective_flux_reconstructor_db));
            
            break;
        }
        case CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS5_JS_HLLC_HLL:
        {
            conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWCNS5_JS_HLLC_HLL(
                "d_convective_flux_reconstructor",
                dim,
                grid_geometry,
                flow_model->getNumberOfEquations(),
                num_species,
                flow_model,
                d_convective_flux_reconstructor_db));
            
            break;
        }
        case CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS5_Z_HLLC_HLL:
        {
            conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWCNS5_Z_HLLC_HLL(
                "d_convective_flux_reconstructor",
                dim,
                grid_geometry,
                flow_model->getNumberOfEquations(),
                num_species,
                flow_model,
                d_convective_flux_reconstructor_db));
            
            break;
        }
        case CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS6_CU_M2_HLLC_HLL:
        {
            conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWCNS6_CU_M2_HLLC_HLL(
                "d_convective_flux_reconstructor",
                dim,
                grid_geometry,
                flow_model->getNumberOfEquations(),
                num_species,
                flow_model,
                d_convective_flux_reconstructor_db));
            
            break;
        }
        case CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS6_LD_HLLC_HLL:
        {
            conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWCNS6_LD_HLLC_HLL(
                "d_convective_flux_reconstructor",
                dim,
                grid_geometry,
                flow_model->getNumberOfEquations(),
                num_species,
                flow_model,
                d_convective_flux_reconstructor_db));
            
            break;
        }
        case CONVECTIVE_FLUX_RECONSTRUCTOR::WCNS6_TEST:
        {
            conv_flux_reconstructor.reset(new ConvectiveFluxReconstructorWCNS6_Test(
                "d_convective_flux_reconstructor",
                dim,
                grid_geometry,
                flow_model->getNumberOfEquations(),
                num_species,
                flow_model,
                d_convective_flux_reconstructor_db));
            
            break;
        }
    }
    
    return conv_flux_reconstructor;
}


/*
 * Print all characteristics of convective flux reconstructor manager.
 */
void
ConvectiveFluxReconstructorManager::printClassData(std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorManager object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorManager: this = "
       << (ConvectiveFluxReconstructorManager *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_convective_flux_reconstructor_type = "
       << d_convective_flux_reconstructor_type
       << std::endl;
}
