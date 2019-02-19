#include "util/mixing_rules/equations_of_state/Noble_Abel_stiffened_gas/EquationOfStateNobleAbelStiffenedGas.hpp"

/*
 * Print all characteristics of the equation of state class.
 */
void
EquationOfStateNobleAbelStiffenedGas::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfStateNobleAbelStiffenedGas object..."
       << std::endl;
       
    os << std::endl;
    
    os << "EquationOfStateNobleAbelStiffenedGas: this = "
       << (EquationOfStateNobleAbelStiffenedGas *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Compute the pressure.
 */
double
EquationOfStateNobleAbelStiffenedGas::getPressure(
    const double* const density,
    const double* const internal_energy,
    const std::vector<const double*>& thermo_properties) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the pressure.
 */
void
EquationOfStateNobleAbelStiffenedGas::computePressure(
    boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the pressure.
 */
void
EquationOfStateNobleAbelStiffenedGas::computePressure(
    boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the pressure.
 */
void
EquationOfStateNobleAbelStiffenedGas::computePressure(
    boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the pressure.
 */
void
EquationOfStateNobleAbelStiffenedGas::computePressure(
    boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the sound speed.
 */
double
EquationOfStateNobleAbelStiffenedGas::getSoundSpeed(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the sound speed.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeSoundSpeed(
    boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the sound speed.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeSoundSpeed(
    boost::shared_ptr<pdat::SideData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the sound speed.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeSoundSpeed(
    boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the sound speed.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeSoundSpeed(
    boost::shared_ptr<pdat::SideData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy.
 */
double
EquationOfStateNobleAbelStiffenedGas::getInternalEnergy(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the specific internal energy.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeInternalEnergy(
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeInternalEnergy(
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeInternalEnergy(
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeInternalEnergy(
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific enthalpy.
 */
double
EquationOfStateNobleAbelStiffenedGas::getEnthalpy(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeEnthalpy(
    boost::shared_ptr<pdat::CellData<double> >& data_enthalpy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeEnthalpy(
    boost::shared_ptr<pdat::SideData<double> >& data_enthalpy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeEnthalpy(
    boost::shared_ptr<pdat::CellData<double> >& data_enthalpy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeEnthalpy(
    boost::shared_ptr<pdat::SideData<double> >& data_enthalpy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the temperature.
 */
double
EquationOfStateNobleAbelStiffenedGas::getTemperature(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the temperature.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeTemperature(
    boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the temperature.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeTemperature(
    boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the temperature.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeTemperature(
    boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the temperature.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeTemperature(
    boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy from temperature.
 */
double
EquationOfStateNobleAbelStiffenedGas::getInternalEnergyFromTemperature(
    const double* const density,
    const double* const temperature,
    const std::vector<const double*>& thermo_properties) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeInternalEnergyFromTemperature(
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeInternalEnergyFromTemperature(
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeInternalEnergyFromTemperature(
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeInternalEnergyFromTemperature(
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the isochoric specific heat capacity.
 */
double
EquationOfStateNobleAbelStiffenedGas::getIsochoricSpecificHeatCapacity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsochoricSpecificHeatCapacity(
    boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsochoricSpecificHeatCapacity(
    boost::shared_ptr<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsochoricSpecificHeatCapacity(
    boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsochoricSpecificHeatCapacity(
    boost::shared_ptr<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the isobaric specific heat capacity.
 */
double
EquationOfStateNobleAbelStiffenedGas::getIsobaricSpecificHeatCapacity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsobaricSpecificHeatCapacity(
    boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsobaricSpecificHeatCapacity(
    boost::shared_ptr<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsobaricSpecificHeatCapacity(
    boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsobaricSpecificHeatCapacity(
    boost::shared_ptr<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
double
EquationOfStateNobleAbelStiffenedGas::getIsochoricPartialInternalEnergyPartialPressure(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsochoricPartialInternalEnergyPartialPressure(
    boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsochoricPartialInternalEnergyPartialPressure(
    boost::shared_ptr<pdat::SideData<double> >& data_partial_internal_energy_partial_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsochoricPartialInternalEnergyPartialPressure(
    boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsochoricPartialInternalEnergyPartialPressure(
    boost::shared_ptr<pdat::SideData<double> >& data_partial_internal_energy_partial_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
double
EquationOfStateNobleAbelStiffenedGas::getIsobaricPartialInternalEnergyPartialDensity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& thermo_properties) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsobaricPartialInternalEnergyPartialDensity(
    boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsobaricPartialInternalEnergyPartialDensity(
    boost::shared_ptr<pdat::SideData<double> >& data_partial_internal_energy_partial_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsobaricPartialInternalEnergyPartialDensity(
    boost::shared_ptr<pdat::CellData<double> >& data_partial_internal_energy_partial_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsobaricPartialInternalEnergyPartialDensity(
    boost::shared_ptr<pdat::SideData<double> >& data_partial_internal_energy_partial_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the density.
 */
double
EquationOfStateNobleAbelStiffenedGas::getDensity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& thermo_properties) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the density.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeDensity(
    boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const std::vector<const double*>& thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the density.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeDensity(
    boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const std::vector<const double*>& thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the density.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeDensity(
    boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_thermo_properties,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the density.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeDensity(
    boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_thermo_properties,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the pressure.
 */
void
EquationOfStateNobleAbelStiffenedGas::computePressure(
    double* const p,
    const double* const rho,
    const double* const epsilon,
    const double& gamma,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_internal_energy,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the pressure.
 */
void
EquationOfStateNobleAbelStiffenedGas::computePressure(
    double* const p,
    const double* const rho,
    const double* const epsilon,
    const double* const gamma,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_internal_energy,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& ghostcell_dims_thermo_properties,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the sound speed.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeSoundSpeed(
    double* const c,
    const double* const rho,
    const double* const p,
    const double& gamma,
    const hier::IntVector& num_ghosts_sound_speed,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& ghostcell_dims_sound_speed,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the sound speed.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeSoundSpeed(
    double* const c,
    const double* const rho,
    const double* const p,
    const double* const gamma,
    const hier::IntVector& num_ghosts_sound_speed,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_sound_speed,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_thermo_properties,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the specific internal energy.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeInternalEnergy(
    double* const epsilon,
    const double* const rho,
    const double* const p,
    const double& gamma,
    const hier::IntVector& num_ghosts_internal_energy,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the specific internal energy.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeInternalEnergy(
    double* const epsilon,
    const double* const rho,
    const double* const p,
    const double* const gamma,
    const hier::IntVector& num_ghosts_internal_energy,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_thermo_properties,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeEnthalpy(
    double* const h,
    const double* const rho,
    const double* const p,
    const double& gamma,
    const hier::IntVector& num_ghosts_enthalpy,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& ghostcell_dims_enthalpy,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the specific enthalpy.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeEnthalpy(
    double* const h,
    const double* const rho,
    const double* const p,
    const double* const gamma,
    const hier::IntVector& num_ghosts_enthalpy,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_enthalpy,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_thermo_properties,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the temperature.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeTemperature(
    double* const T,
    const double* const rho,
    const double* const p,
    const double& gamma,
    const double& c_v,
    const hier::IntVector& num_ghosts_temperature,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& ghostcell_dims_temperature,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the temperature.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeTemperature(
    double* const T,
    const double* const rho,
    const double* const p,
    const double* const gamma,
    const double* const c_v,
    const hier::IntVector& num_ghosts_temperature,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_temperature,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_thermo_properties,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeInternalEnergyFromTemperature(
    double* const epsilon,
    const double* const T,
    const double& c_v,
    const hier::IntVector& num_ghosts_internal_energy,
    const hier::IntVector& num_ghosts_temperature,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& ghostcell_dims_temperature,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the specific internal energy from temperature.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeInternalEnergyFromTemperature(
    double* const epsilon,
    const double* const T,
    const double* const c_v,
    const hier::IntVector& num_ghosts_internal_energy,
    const hier::IntVector& num_ghosts_temperature,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_internal_energy,
    const hier::IntVector& ghostcell_dims_temperature,
    const hier::IntVector& ghostcell_dims_thermo_properties,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsochoricSpecificHeatCapacity(
    double* const c_v,
    const double& c_v_src,
    const hier::IntVector& num_ghosts_isochoric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the isochoric specific heat capacity.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsochoricSpecificHeatCapacity(
    double* const c_v,
    const double* const c_v_src,
    const hier::IntVector& num_ghosts_isochoric_specific_heat_capacity,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_isochoric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_thermo_properties,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsobaricSpecificHeatCapacity(
    double* const c_p,
    const double& c_p_src,
    const hier::IntVector& num_ghosts_isobaric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the isobaric specific heat capacity.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsobaricSpecificHeatCapacity(
    double* const c_p,
    const double* const c_p_src,
    const hier::IntVector& num_ghosts_isobaric_specific_heat_capacity,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_isobaric_specific_heat_capacity,
    const hier::IntVector& ghostcell_dims_thermo_properties,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsochoricPartialInternalEnergyPartialPressure(
    double* const xi,
    const double& gamma,
    const hier::IntVector& num_ghosts_partial_internal_energy_partial_pressure,
    const hier::IntVector& ghostcell_dims_partial_internal_energy_partial_pressure,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the partial derivative of internal energy w.r.t. pressure under constant density.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsochoricPartialInternalEnergyPartialPressure(
    double* const xi,
    const double* const gamma,
    const hier::IntVector& num_ghosts_partial_internal_energy_partial_pressure,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_partial_internal_energy_partial_pressure,
    const hier::IntVector& ghostcell_dims_thermo_properties,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the partial derivative of internal energy w.r.t. density under constant pressure.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeIsobaricPartialInternalEnergyPartialDensity(
    double* const delta,
    const hier::IntVector& num_ghosts_partial_internal_energy_partial_density,
    const hier::IntVector& ghostcell_dims_partial_internal_energy_partial_density,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the density.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeDensity(
    double* const rho,
    const double* const p,
    const double* const T,
    const double& R,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_temperature,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_temperature,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}


// Need to redefine this helper function.
/*
 * Compute the density.
 */
void
EquationOfStateNobleAbelStiffenedGas::computeDensity(
    double* const rho,
    const double* const p,
    const double* const T,
    const double* const R,
    const hier::IntVector& num_ghosts_density,
    const hier::IntVector& num_ghosts_pressure,
    const hier::IntVector& num_ghosts_temperature,
    const hier::IntVector& num_ghosts_thermo_properties,
    const hier::IntVector& ghostcell_dims_density,
    const hier::IntVector& ghostcell_dims_pressure,
    const hier::IntVector& ghostcell_dims_temperature,
    const hier::IntVector& ghostcell_dims_thermo_properties,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims) const
{
    // NEED IMPLEMENTATION!
}
