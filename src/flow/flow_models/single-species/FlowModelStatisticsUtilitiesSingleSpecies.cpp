#include "flow/flow_models/single-species/FlowModelStatisticsUtilitiesSingleSpecies.hpp"

#include "SAMRAI/hier/FlattenedHierarchy.h"

#include <fstream>

/*
 * Output names of statistical quantities to output to a file.
 */
void
FlowModelStatisticsUtilitiesSingleSpecies::outputStatisticalQuantitiesNames(
    const std::string& stat_dump_filename)
{
    NULL_USE(stat_dump_filename);
}


/*
 * Output statisitcal quantities to a file.
 */
void
FlowModelStatisticsUtilitiesSingleSpecies::outputStatisticalQuantities(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double output_time)
{
    NULL_USE(stat_dump_filename);
    NULL_USE(patch_hierarchy);
    NULL_USE(data_context);
    NULL_USE(output_time);
}