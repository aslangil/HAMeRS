#ifndef FLOW_MODEL_FOUR_EQN_CONSERVATIVE_HPP
#define FLOW_MODEL_FOUR_EQN_CONSERVATIVE_HPP

#include "flow/flow_models/FlowModel.hpp"
#include "util/mixing_rules/equations_of_mass_diffusivity/EquationOfMassDiffusivityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_shear_viscosity/EquationOfShearViscosityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_bulk_viscosity/EquationOfBulkViscosityMixingRulesManager.hpp"
#include "util/mixing_rules/equations_of_thermal_conductivity/EquationOfThermalConductivityMixingRulesManager.hpp"

class FlowModelFourEqnConservative: public FlowModel
{
    public:
        FlowModelFourEqnConservative(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<tbox::Database>& flow_model_db);
        
        ~FlowModelFourEqnConservative() {}
        
        /*
         * Print all characteristics of the flow model class.
         */
        void printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the flow model class into the restart database.
         */
        void putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Register the conservative variables.
         */
        void
        registerConservativeVariables(
            RungeKuttaLevelIntegrator* integrator,
            const hier::IntVector& num_ghosts,
            const hier::IntVector& num_ghosts_intermediate);
        
        /*
         * Get the names of conservative variables.
         */
        std::vector<std::string> getNamesOfConservativeVariables(bool have_underscores = false);
        
        /*
         * Get the names of primitive variables.
         */
        std::vector<std::string> getNamesOfPrimitiveVariables(bool have_underscores = false);
        
        /*
         * Get the variable types of conservative variables.
         */
        std::vector<std::string> getVariableTypesOfConservativeVariables();
        
        /*
         * Get the variable types of primitive variables.
         */
        std::vector<std::string> getVariableTypesOfPrimitiveVariables();
        
        /*
         * Get the conservative variables.
         */
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > >
        getConservativeVariables();
        
        /*
         * Register a patch with a data context.
         */
        void
        registerPatchWithDataContext(
            const hier::Patch& patch,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Register different derived variables in the registered patch. The derived variables to be registered
         * are given as entires in a map of the variable name to the number of sub-ghost cells required.
         * If the variable to be registered is one of the conservative variable, the corresponding entry
         * in the map is ignored.
         */
        void
        registerDerivedVariables(
            const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data);
        
        /*
         * Unregister the registered patch. The registered data context and the cell data of all derived variables in
         * the patch are cleared.
         */
        void unregisterPatch();
        
        /*
         * Allocate memory for cell data of different registered derived variables.
         */
        void allocateMemoryForDerivedCellData();
        
        /*
         * Compute the cell data of different registered derived variables with the registered data context.
         */
        void
        computeDerivedCellData();
        
        /*
         * Get the cell data of one cell variable in the registered patch.
         */
        boost::shared_ptr<pdat::CellData<double> >
        getCellData(const std::string& variable_key);
        
        /*
         * Get the cell data of different cell variables in the registered patch.
         */
        std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getCellData(const std::vector<std::string>& variable_keys);
        
        /*
         * Get the cell data of species cell variables in the registered patch.
         */
        std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getSpeciesCellData(const std::string& variable_key);
        
        /*
         * Fill the cell data of conservative variables in the interior box with value zero.
         */
        void
        fillCellDataOfConservativeVariablesWithZero();
        
        /*
         * Update the cell data of conservative variables in the interior box after time advancement.
         */
        void
        updateCellDataOfConservativeVariables();
        
        /*
         * Get the cell data of the conservative variables in the registered patch.
         */
        std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getCellDataOfConservativeVariables();
        
        /*
         * Get the cell data of the primitive variables in the registered patch.
         */
        std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getCellDataOfPrimitiveVariables();
        
        /*
         * Compute derived plot quantities registered with the VisIt data writers from data that
         * is maintained on each patch in the hierarchy.
         */
        bool
        packDerivedDataIntoDoubleBuffer(
            double* buffer,
            const hier::Patch& patch,
            const hier::Box& region,
            const std::string& variable_name,
            int depth_id,
            double simulation_time) const;
        
        /*
         * Register the plotting quantities.
         */
#ifdef HAVE_HDF5
        void
        registerPlotQuantities(
            const boost::shared_ptr<ExtendedVisItDataWriter>& visit_writer);
#endif
        
    private:
        /*
         * Set the number of sub-ghost cells of a variable.
         * This function can be called recursively if the variables are computed recursively.
         */
        void
        setNumberOfSubGhosts(
            const hier::IntVector& num_subghosts,
            const std::string& variable_name,
            const std::string& parent_variable_name);
        
        /*
         * Set the ghost boxes of derived cell variables.
         */
        void
        setDerivedCellVariableGhostBoxes();
        
        /*
         * Get the cell data of partial densities in the registered patch.
         */
        boost::shared_ptr<pdat::CellData<double> >
        getCellDataOfPartialDensities();
        
        /*
         * Get the cell data of momentum in the registered patch.
         */
        boost::shared_ptr<pdat::CellData<double> >
        getCellDataOfMomentum();
        
        /*
         * Get the cell data of total energy in the registered patch.
         */
        boost::shared_ptr<pdat::CellData<double> >
        getCellDataOfTotalEnergy();
        
        /*
         * Compute the cell data of density in the registered patch.
         */
        void computeCellDataOfDensity(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of mass fractions with density in the registered patch.
         */
        void computeCellDataOfMassFractionsWithDensity(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of mole fractions with mass fractions in the registered patch.
         */
        void computeCellDataOfMoleFractionsWithMassFractions(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of velocity with density in the registered patch.
         */
        void computeCellDataOfVelocityWithDensity(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of internal energy with density and velocity in the registered
         * patch.
         */
        void computeCellDataOfInternalEnergyWithDensityAndVelocity(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of pressure with density, mass fractions and internal energy in the registered patch.
         */
        void computeCellDataOfPressureWithDensityMassFractionsAndInternalEnergy(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of sound speed with density, mass fractions and pressure in the registered patch.
         */
        void computeCellDataOfSoundSpeedWithDensityMassFractionsAndPressure(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of temperature with density, mass fractions and pressure in the registered patch.
         */
        void computeCellDataOfTemperatureWithDensityMassFractionsAndPressure(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of convective flux with velocity and pressure in the registered patch.
         */
        void computeCellDataOfConvectiveFluxWithVelocityAndPressure(
            const DIRECTION::TYPE& direction,
            const hier::Box& domain);
        
        /*
         * Compute the cell data of maximum wave speed with velocity and sound speed in the registered patch.
         */
        void computeCellDataOfMaxWaveSpeedWithVelocityAndSoundSpeed(
            const DIRECTION::TYPE& direction,
            const hier::Box& domain);
        
        /*
         * Compute the cell data of maximum diffusivity with density, mass fractions, pressure and temperature in the
         * registered patch.
         */
        void computeCellDataOfMaxDiffusivityWithDensityMassFractionsPressureAndTemperature(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of species densities with pressure and temperature in the registered patch.
         */
        void computeCellDataOfSpeciesDensitiesWithPressureAndTemperature(
            const hier::Box& domain);
        
        /*
         * Compute the cell data of species enthalpies with species densities and pressure in the registered patch.
         */
        void computeCellDataOfSpeciesEnthalpiesWithSpeciesDensitiesAndPressure(
            const hier::Box& domain);
        
        /*
         * boost::shared_ptr to registered conservative variables.
         */
        static boost::shared_ptr<pdat::CellVariable<double> > s_variable_partial_densities;
        static boost::shared_ptr<pdat::CellVariable<double> > s_variable_momentum;
        static boost::shared_ptr<pdat::CellVariable<double> > s_variable_total_energy;
        
        /*
         * Number of sub-ghost cells of derived cell data.
         */
        hier::IntVector d_num_subghosts_density;
        hier::IntVector d_num_subghosts_mass_fractions;
        hier::IntVector d_num_subghosts_mole_fractions;
        hier::IntVector d_num_subghosts_velocity;
        hier::IntVector d_num_subghosts_internal_energy;
        hier::IntVector d_num_subghosts_pressure;
        hier::IntVector d_num_subghosts_sound_speed;
        hier::IntVector d_num_subghosts_temperature;
        hier::IntVector d_num_subghosts_convective_flux_x;
        hier::IntVector d_num_subghosts_convective_flux_y;
        hier::IntVector d_num_subghosts_convective_flux_z;
        hier::IntVector d_num_subghosts_max_wave_speed_x;
        hier::IntVector d_num_subghosts_max_wave_speed_y;
        hier::IntVector d_num_subghosts_max_wave_speed_z;
        hier::IntVector d_num_subghosts_max_diffusivity;
        hier::IntVector d_num_subghosts_species_densities;
        hier::IntVector d_num_subghosts_species_enthalpies;
        
        /*
         * Boxes with sub-ghost cells of derived cell data.
         */
        hier::Box d_subghost_box_density;
        hier::Box d_subghost_box_mass_fractions;
        hier::Box d_subghost_box_mole_fractions;
        hier::Box d_subghost_box_velocity;
        hier::Box d_subghost_box_internal_energy;
        hier::Box d_subghost_box_pressure;
        hier::Box d_subghost_box_sound_speed;
        hier::Box d_subghost_box_temperature;
        hier::Box d_subghost_box_convective_flux_x;
        hier::Box d_subghost_box_convective_flux_y;
        hier::Box d_subghost_box_convective_flux_z;
        hier::Box d_subghost_box_max_wave_speed_x;
        hier::Box d_subghost_box_max_wave_speed_y;
        hier::Box d_subghost_box_max_wave_speed_z;
        hier::Box d_subghost_box_max_diffusivity;
        hier::Box d_subghost_box_species_densities;
        hier::Box d_subghost_box_species_enthalpies;
        
        /*
         * Dimensions of boxes with sub-ghost cells of derived cell data.
         */
        hier::IntVector d_subghostcell_dims_density;
        hier::IntVector d_subghostcell_dims_mass_fractions;
        hier::IntVector d_subghostcell_dims_mole_fractions;
        hier::IntVector d_subghostcell_dims_velocity;
        hier::IntVector d_subghostcell_dims_internal_energy;
        hier::IntVector d_subghostcell_dims_pressure;
        hier::IntVector d_subghostcell_dims_sound_speed;
        hier::IntVector d_subghostcell_dims_temperature;
        hier::IntVector d_subghostcell_dims_convective_flux_x;
        hier::IntVector d_subghostcell_dims_convective_flux_y;
        hier::IntVector d_subghostcell_dims_convective_flux_z;
        hier::IntVector d_subghostcell_dims_max_wave_speed_x;
        hier::IntVector d_subghostcell_dims_max_wave_speed_y;
        hier::IntVector d_subghostcell_dims_max_wave_speed_z;
        hier::IntVector d_subghostcell_dims_max_diffusivity;
        hier::IntVector d_subghostcell_dims_species_densities;
        hier::IntVector d_subghostcell_dims_species_enthalpies;
        
        /*
         * boost::shared_ptr to derived cell data.
         */
        boost::shared_ptr<pdat::CellData<double> > d_data_density;
        boost::shared_ptr<pdat::CellData<double> > d_data_mass_fractions;
        boost::shared_ptr<pdat::CellData<double> > d_data_mole_fractions;
        boost::shared_ptr<pdat::CellData<double> > d_data_velocity;
        boost::shared_ptr<pdat::CellData<double> > d_data_internal_energy;
        boost::shared_ptr<pdat::CellData<double> > d_data_pressure;
        boost::shared_ptr<pdat::CellData<double> > d_data_sound_speed;
        boost::shared_ptr<pdat::CellData<double> > d_data_temperature;
        boost::shared_ptr<pdat::CellData<double> > d_data_convective_flux_x;
        boost::shared_ptr<pdat::CellData<double> > d_data_convective_flux_y;
        boost::shared_ptr<pdat::CellData<double> > d_data_convective_flux_z;
        boost::shared_ptr<pdat::CellData<double> > d_data_max_wave_speed_x;
        boost::shared_ptr<pdat::CellData<double> > d_data_max_wave_speed_y;
        boost::shared_ptr<pdat::CellData<double> > d_data_max_wave_speed_z;
        boost::shared_ptr<pdat::CellData<double> > d_data_max_diffusivity;
        std::vector<boost::shared_ptr<pdat::CellData<double> > > d_data_species_densities;
        std::vector<boost::shared_ptr<pdat::CellData<double> > > d_data_species_enthalpies;
        
        /*
         * Whether derived cell data is computed.
         */
        bool d_cell_data_computed_density;
        bool d_cell_data_computed_mass_fractions;
        bool d_cell_data_computed_mole_fractions;
        bool d_cell_data_computed_velocity;
        bool d_cell_data_computed_internal_energy;
        bool d_cell_data_computed_pressure;
        bool d_cell_data_computed_sound_speed;
        bool d_cell_data_computed_temperature;
        bool d_cell_data_computed_convective_flux_x;
        bool d_cell_data_computed_convective_flux_y;
        bool d_cell_data_computed_convective_flux_z;
        bool d_cell_data_computed_max_wave_speed_x;
        bool d_cell_data_computed_max_wave_speed_y;
        bool d_cell_data_computed_max_wave_speed_z;
        bool d_cell_data_computed_max_diffusivity;
        bool d_cell_data_computed_species_densities;
        bool d_cell_data_computed_species_enthalpies;
        
        /*
         * A string variable to describe the equation of mass diffusivity used.
         */
        std::string d_equation_of_mass_diffusivity_str;
        
        /*
         * A string variable to describe the equation of shear viscosity used.
         */
        std::string d_equation_of_shear_viscosity_str;
        
        /*
         * A string variable to describe the equation of bulk viscosity used.
         */
        std::string d_equation_of_bulk_viscosity_str;
        
        /*
         * A string variable to describe the equation of thermal conductivity used.
         */
        std::string d_equation_of_thermal_conductivity_str;
        
        /*
         * boost::shared_ptr to EquationOfMassDiffusivityMixingRules.
         */
        boost::shared_ptr<EquationOfMassDiffusivityMixingRules>
            d_equation_of_mass_diffusivity_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfMassDiffusivityMixingRulesManager.
         */
        boost::shared_ptr<EquationOfMassDiffusivityMixingRulesManager>
            d_equation_of_mass_diffusivity_mixing_rules_manager;
        
        /*
         * boost::shared_ptr to EquationOfShearViscosityMixingRules.
         */
        boost::shared_ptr<EquationOfShearViscosityMixingRules>
            d_equation_of_shear_viscosity_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfShearViscosityMixingRulesManager.
         */
        boost::shared_ptr<EquationOfShearViscosityMixingRulesManager>
            d_equation_of_shear_viscosity_mixing_rules_manager;
        
        /*
         * boost::shared_ptr to EquationOfBulkViscosityMixingRules.
         */
        boost::shared_ptr<EquationOfBulkViscosityMixingRules>
            d_equation_of_bulk_viscosity_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfBulkViscosityMixingRulesManager.
         */
        boost::shared_ptr<EquationOfBulkViscosityMixingRulesManager>
            d_equation_of_bulk_viscosity_mixing_rules_manager;
        
        /*
         * boost::shared_ptr to EquationOfThermalConductivityMixingRules.
         */
        boost::shared_ptr<EquationOfThermalConductivityMixingRules>
            d_equation_of_thermal_conductivity_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfThermalConductivityMixingRulesManager.
         */
        boost::shared_ptr<EquationOfThermalConductivityMixingRulesManager>
            d_equation_of_thermal_conductivity_mixing_rules_manager;
        
};

#endif /* FLOW_MODEL_FOUR_EQN_CONSERVATIVE_HPP */
