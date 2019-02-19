#include "util/mixing_rules/equations_of_state/Noble_Abel_stiffened_gas/EquationOfStateMixingRulesNobleAbelStiffenedGas.hpp"

EquationOfStateMixingRulesNobleAbelStiffenedGas::EquationOfStateMixingRulesNobleAbelStiffenedGas(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const MIXING_CLOSURE_MODEL::TYPE& mixing_closure_model,
    const boost::shared_ptr<tbox::Database>& equation_of_state_mixing_rules_db):
        EquationOfStateMixingRules(
            object_name,
            dim,
            num_species,
            mixing_closure_model,
            equation_of_state_mixing_rules_db)
{
    d_equation_of_state.reset(new EquationOfStateNobleAbelStiffenedGas(
        "d_equation_of_state",
        dim));
    
    /*
     * Get the ratio of specific heats of each species from the database.
     */
    
    if (equation_of_state_mixing_rules_db->keyExists("species_gamma"))
    {
        size_t species_gamma_array_size = equation_of_state_mixing_rules_db->getArraySize("species_gamma");
        if (static_cast<int>(species_gamma_array_size) == d_num_species)
        {
            d_species_gamma = equation_of_state_mixing_rules_db->getDoubleVector("species_gamma");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_gamma' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_gamma"))
    {
        size_t species_gamma_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_gamma");
        if (static_cast<int>(species_gamma_array_size) == d_num_species)
        {
            d_species_gamma = equation_of_state_mixing_rules_db->getDoubleVector("d_species_gamma");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_gamma' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_gamma'/'d_species_gamma'"
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    /*
     * Get the Noble-Abel stiffened gas coefficients of each species from the database.
     */
    
    if (equation_of_state_mixing_rules_db->keyExists("species_p_inf"))
    {
        size_t species_p_inf_array_size = equation_of_state_mixing_rules_db->getArraySize("species_p_inf");
        if (static_cast<int>(species_p_inf_array_size) == d_num_species)
        {
            d_species_p_inf = equation_of_state_mixing_rules_db->getDoubleVector("species_p_inf");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_p_inf' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_p_inf"))
    {
        size_t species_p_inf_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_p_inf");
        if (static_cast<int>(species_p_inf_array_size) == d_num_species)
        {
            d_species_p_inf = equation_of_state_mixing_rules_db->getDoubleVector("d_species_p_inf");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_p_inf' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_p_inf'/'d_species_p_inf'"
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    if (equation_of_state_mixing_rules_db->keyExists("species_q"))
    {
        size_t species_q_array_size = equation_of_state_mixing_rules_db->getArraySize("species_q");
        if (static_cast<int>(species_q_array_size) == d_num_species)
        {
            d_species_q = equation_of_state_mixing_rules_db->getDoubleVector("species_q");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_q' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_q"))
    {
        size_t species_q_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_q");
        if (static_cast<int>(species_q_array_size) == d_num_species)
        {
            d_species_q = equation_of_state_mixing_rules_db->getDoubleVector("d_species_q");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_q' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_q'/'d_species_q'"
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    if (equation_of_state_mixing_rules_db->keyExists("species_q_prime"))
    {
        size_t species_q_prime_array_size = equation_of_state_mixing_rules_db->getArraySize("species_q_prime");
        if (static_cast<int>(species_q_prime_array_size) == d_num_species)
        {
            d_species_q_prime = equation_of_state_mixing_rules_db->getDoubleVector("species_q_prime");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_q_prime' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_q_prime"))
    {
        size_t species_q_prime_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_q_prime");
        if (static_cast<int>(species_q_prime_array_size) == d_num_species)
        {
            d_species_q_prime = equation_of_state_mixing_rules_db->getDoubleVector("d_species_q_prime");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_q_prime' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_q_prime'/'d_species_q_prime'"
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    if (equation_of_state_mixing_rules_db->keyExists("species_b"))
    {
        size_t species_b_array_size = equation_of_state_mixing_rules_db->getArraySize("species_b");
        if (static_cast<int>(species_b_array_size) == d_num_species)
        {
            d_species_b = equation_of_state_mixing_rules_db->getDoubleVector("species_b");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_b' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_b"))
    {
        size_t species_b_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_b");
        if (static_cast<int>(species_b_array_size) == d_num_species)
        {
            d_species_b = equation_of_state_mixing_rules_db->getDoubleVector("d_species_b");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_b' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_b'/'d_species_b'"
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    /*
     * Get the specific heat at constant volume of each species from the database.
     */
    
    if (equation_of_state_mixing_rules_db->keyExists("species_c_v"))
    {
        size_t species_c_v_array_size = equation_of_state_mixing_rules_db->getArraySize("species_c_v");
        if (static_cast<int>(species_c_v_array_size) == d_num_species)
        {
            d_species_c_v = equation_of_state_mixing_rules_db->getDoubleVector("species_c_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_c_v' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_c_v"))
    {
        size_t species_c_v_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_c_v");
        if (static_cast<int>(species_c_v_array_size) == d_num_species)
        {
            d_species_c_v = equation_of_state_mixing_rules_db->getDoubleVector("d_species_c_v");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_c_v' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_c_v'/'d_species_c_v'"
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    /*
     * Get the molcular weight of each species from the database.
     */
    
    if (equation_of_state_mixing_rules_db->keyExists("species_M"))
    {
        size_t species_M_array_size = equation_of_state_mixing_rules_db->getArraySize("species_M");
        if (static_cast<int>(species_M_array_size) == d_num_species)
        {
            d_species_M = equation_of_state_mixing_rules_db->getDoubleVector("species_M");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'species_M' entries must be equal to 'num_species'."
                << std::endl);
        }
    }
    else if (equation_of_state_mixing_rules_db->keyExists("d_species_M"))
    {
        size_t species_M_array_size = equation_of_state_mixing_rules_db->getArraySize("d_species_M");
        if (static_cast<int>(species_M_array_size) == d_num_species)
        {
            d_species_M = equation_of_state_mixing_rules_db->getDoubleVector("d_species_M");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "number of 'd_species_M' entries must be equal to 'd_num_species'."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Key data 'species_M'/'d_species_M'"
            << "not found in data for Equation_of_state_mixing_rules."
            << std::endl);
    }
    
    /*
     * Compute the specific heat at constant pressure of each species.
     */
    
    d_species_c_p.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        d_species_c_p.push_back(d_species_gamma[si]*d_species_c_v[si]);
    }
}


/*
 * Print all characteristics of the equation of state class.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfStateMixingRulesNobleAbelStiffenedGas object..."
       << std::endl;
    
    os << std::endl;
    os << "EquationOfStateMixingRulesNobleAbelStiffenedGas: this = "
       << (EquationOfStateMixingRulesNobleAbelStiffenedGas *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    
    os << "d_mixing_closure_model = "
       << d_mixing_closure_model
       << std::endl;
    
    /*
     * Print the ratio of specific heats of each species.
     */
    
    os << "d_species_gamma = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_gamma[si] << ", ";
    }
    os << d_species_gamma[d_num_species - 1];
    os << std::endl;
    
    /*
     * Print the Noble-Abel stiffened gas coefficients of each species.
     */
    
    os << "d_species_p_inf = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_p_inf[si] << ", ";
    }
    os << d_species_p_inf[d_num_species - 1];
    os << std::endl;
    
    os << "d_species_q = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_q[si] << ", ";
    }
    os << d_species_q[d_num_species - 1];
    os << std::endl;
    
    os << "d_species_q_prime = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_q_prime[si] << ", ";
    }
    os << d_species_q_prime[d_num_species - 1];
    os << std::endl;
    
    os << "d_species_b = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_b[si] << ", ";
    }
    os << d_species_b[d_num_species - 1];
    os << std::endl;
    
    /*
     * Print the specific heat at constant volume of each species.
     */
    
    os << "d_species_c_v = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_c_v[si] << ", ";
    }
    os << d_species_c_v[d_num_species - 1];
    os << std::endl;
    
    /*
     * Print the molecular weight of each species.
     */
    
    os << "d_species_M = ";
    for (int si = 0; si < d_num_species - 1; si++)
    {
        os << d_species_M[si] << ", ";
    }
    os << d_species_M[d_num_species - 1];
    os << std::endl;
}


/*
 * Put the characteristics of the equation of state mixing rules class into the restart
 * database.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putDoubleVector("d_species_gamma", d_species_gamma);
    restart_db->putDoubleVector("d_species_p_inf", d_species_p_inf);
    restart_db->putDoubleVector("d_species_q", d_species_q);
    restart_db->putDoubleVector("d_species_q_prime", d_species_q_prime);
    restart_db->putDoubleVector("d_species_b", d_species_b);
    restart_db->putDoubleVector("d_species_c_v", d_species_c_v);
    restart_db->putDoubleVector("d_species_M", d_species_M);
}


/*
 * Compute the pressure of the mixture with isothermal and isobaric equilibria assumptions.
 */
double
EquationOfStateMixingRulesNobleAbelStiffenedGas::getPressure(
    const double* const density,
    const double* const internal_energy,
    const std::vector<const double*>& mass_fractions) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the pressure of the mixture with isothermal and isobaric equilibria assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computePressure(
    boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the pressure of the mixture with isothermal and isobaric equilibria assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computePressure(
    boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the pressure of the mixture with isobaric equilibrium assumption.
 */
double
EquationOfStateMixingRulesNobleAbelStiffenedGas::getPressure(
    const double* const density,
    const double* const internal_energy,
    const std::vector<const double*>& mass_fractions,
    const std::vector<const double*>& volume_fractions) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the pressure of the mixture with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computePressure(
    boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the pressure of the mixture with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computePressure(
    boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the sound speed of the mixture with isothermal and isobaric equilibria assumptions.
 */
double
EquationOfStateMixingRulesNobleAbelStiffenedGas::getSoundSpeed(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the sound speed of the mixture with isothermal and isobaric equilibria assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeSoundSpeed(
    boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the sound speed of the mixture with isothermal and isobaric equilibria assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeSoundSpeed(
    boost::shared_ptr<pdat::SideData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the sound speed of the mixture with isobaric equilibrium assumption.
 */
double
EquationOfStateMixingRulesNobleAbelStiffenedGas::getSoundSpeed(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions,
    const std::vector<const double*>& volume_fractions) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the sound speed of the mixture with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeSoundSpeed(
    boost::shared_ptr<pdat::CellData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the sound speed of the mixture with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeSoundSpeed(
    boost::shared_ptr<pdat::SideData<double> >& data_sound_speed,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy of the mixture with isothermal and isobaric equilibria assumptions.
 */
double
EquationOfStateMixingRulesNobleAbelStiffenedGas::getInternalEnergy(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the specific internal energy of the mixture with isothermal and isobaric equilibria assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeInternalEnergy(
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy of the mixture with isothermal and isobaric equilibria assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeInternalEnergy(
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
 */
double
EquationOfStateMixingRulesNobleAbelStiffenedGas::getInternalEnergy(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions,
    const std::vector<const double*>& volume_fractions) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeInternalEnergy(
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::CellData<double> >& data_volume_fractions,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy of the mixture with isobaric equilibrium assumption.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeInternalEnergy(
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    const boost::shared_ptr<pdat::SideData<double> >& data_volume_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the temperature of the mixture with isothermal and isobaric equilibria assumptions.
 */
double
EquationOfStateMixingRulesNobleAbelStiffenedGas::getTemperature(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the temperature of the mixture with isothermal and isobaric equilibria assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeTemperature(
    boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the temperature of the mixture with isothermal and isobaric equilibria assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeTemperature(
    boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy of the mixture from temperature with isothermal
 * and isobaric equilibria assumptions.
 */
double
EquationOfStateMixingRulesNobleAbelStiffenedGas::getInternalEnergyFromTemperature(
    const double* const density,
    const double* const temperature,
    const std::vector<const double*>& mass_fractions) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the specific internal energy of the mixture from temperature with isothermal
 * and isobaric equilibria assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeInternalEnergyFromTemperature(
    boost::shared_ptr<pdat::CellData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the specific internal energy of the mixture from temperature with isothermal
 * and isobaric equilibria assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeInternalEnergyFromTemperature(
    boost::shared_ptr<pdat::SideData<double> >& data_internal_energy,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric equilibria
 * assumptions.
 */
double
EquationOfStateMixingRulesNobleAbelStiffenedGas::getIsochoricSpecificHeatCapacity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric equilibria
 * assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeIsochoricSpecificHeatCapacity(
    boost::shared_ptr<pdat::CellData<double> >& data_isochoric_specific_heat_capacity,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the isochoric specific heat capacity of mixture with isothermal and isobaric equilibria
 * assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeIsochoricSpecificHeatCapacity(
    boost::shared_ptr<pdat::SideData<double> >& data_isochoric_specific_heat_capacity,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric equilibria
 * assumptions.
 */
double
EquationOfStateMixingRulesNobleAbelStiffenedGas::getIsobaricSpecificHeatCapacity(
    const double* const density,
    const double* const pressure,
    const std::vector<const double*>& mass_fractions) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric equilibria
 * assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeIsobaricSpecificHeatCapacity(
    boost::shared_ptr<pdat::CellData<double> >& data_isobaric_specific_heat_capacity,
    const boost::shared_ptr<pdat::CellData<double> >& data_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the isobaric specific heat capacity of mixture with isothermal and isobaric equilibria
 * assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeIsobaricSpecificHeatCapacity(
    boost::shared_ptr<pdat::SideData<double> >& data_isobaric_specific_heat_capacity,
    const boost::shared_ptr<pdat::SideData<double> >& data_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the density of mixture with isothermal and isobaric equilibria assumptions.
 */
double
EquationOfStateMixingRulesNobleAbelStiffenedGas::getMixtureDensity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& mass_fractions) const
{
    // NEED IMPLEMENTATION!
    return double(0);
}


/*
 * Compute the density of mixture with isothermal and isobaric equilibria assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeMixtureDensity(
    boost::shared_ptr<pdat::CellData<double> >& data_mixture_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Compute the density of mixture with isothermal and isobaric equilibria assumptions.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::computeMixtureDensity(
    boost::shared_ptr<pdat::SideData<double> >& data_mixture_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_pressure,
    const boost::shared_ptr<pdat::SideData<double> >& data_temperature,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
    // NEED IMPLEMENTATION!
}


/*
 * Get the thermodynamic properties of a species.
 */
void
EquationOfStateMixingRulesNobleAbelStiffenedGas::getSpeciesThermodynamicProperties(
    std::vector<double*>& species_thermo_properties,
    const int species_index) const
{
    // NEED IMPLEMENTATION!
}
