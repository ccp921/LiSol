! This module is part of LiSol
! This module is used for reading in mechanism through Cantera,
! and does the following calculation by self, which includes k-
! inetic, thermodynamic, and transport properties.

!!!!! NOTE the units here !!!!!

module lisol_chemistry_mod
    use cantera
    implicit none

    character(len=300), parameter :: mech_path = './mech' 
    character(len=200) :: mech_name
    type(phase_t) gas

    !=========================
    ! variables for elements
    !=========================
    integer(kind=4) :: nelem  !< number of elements
    character(len=2), allocatable, dimension(:) :: elem_name   !< name strings of elements 
    !=========================
    ! variables for species
    !=========================
    integer(kind=4) :: nspec  !< number of species
    character(len=20), allocatable, dimension(:) :: spec_name  !< name strings of species
    integer(kind=4), allocatable, dimension(:,:) :: elem_num   !< element composition of each species, e.g. CO2 has 'C:1, O:2', SIZE: nspec*nelem
    integer(kind=4), parameter :: spec_therm_num = 17          !< length of array spec_therm for single species, 17 (lowT, midT, highT, coeffs in [lowT,midT], and coeffs in [midT,highT])
    real(kind=8), allocatable, dimension(:,:) :: spec_therm    !< thermodynamic data, NASA7 coefficients, SIZE: nspec*spec_therm_num
    integer(kind=4), parameter :: spec_trans_num = 6           !< length of array spec_trans for single species, 6 (geometry(0-single atom, 1-linear, 2-nonlinear), LJ well depth(K), LJ diameter(A), dipole, polarizability, and rotational relaxation)
    real(kind=8), allocatable, dimension(:,:) :: spec_trans    !< transport data: 1 - geometry (0-single atom, 1-linear, 2-nonlinear)
                                                               !<                 2 - LJ well depth [K]
                                                               !<                 3 - LJ diameter [A]
                                                               !<                 4 - dipole [Debye]
                                                               !<                 5 - polarizability [A^3]
                                                               !<                 6 - rotational relaxation [-]
                                                               !< SIZE: nspec*6
    real(kind=8), allocatable, dimension(:,:) :: spec_trans_ISO !< transport data in ISO units, well depth[J], diameter[m], dipole[Coulomb-m], polarizability[m^3]
    character(len=50), parameter :: inter_path = './exchange_data'     !< used to get thermo and transport data from Cantera, through Python scripts output and data read-in
    character(len=50), parameter :: spec_file_suffix = '_spec.dat'     !< element atom number in specific species
    character(len=50), parameter :: thermo_file_suffix = '_thermo.dat' !< thermo data of specific species
    character(len=50), parameter :: trans_file_suffix = '_trans.dat'   !< trans data of specific species
    character(len=50), parameter :: py_file = './py2fort.py'           !< file name of the python to fortran script
    character(len=300), parameter :: fort2py = './fort2py.dat'         !< transfer the mechanism file name
    !=========================
    ! variables for reactions
    !=========================
    integer(kind=4) :: nreac  !< number of reactions
    character(len=200), allocatable, dimension(:) :: reac_name !< strings of reaction equations, SIZE: nreac*length
    integer(kind=4), allocatable, dimension(:,:) :: spec_loc   !< if specific species occurs in the reactants or products of specific reaction, 0-NO, 1-YES, SIZE: nspec*nreac*2
                                                                ! can be obtained from Cantera directly from functions reactantStoichCoeff and productStoichCoeff
    integer(kind=4), allocatable, dimension(:,:) :: reac_type    !< types of each reaction, SIZE: nreac*5, refering to variables defined below
    logical, allocatable, dimension(:) :: isRevers !< If the reactions are reversible
    integer(kind=4), allocatable, dimension(:,:) :: reac_coeffs !< stoichiometric coefficients of each species in the reactants among each reaction, SIZE: nspec*nreac
    integer(kind=4), allocatable, dimension(:,:) :: prod_coeffs !< stoichiometric coefficients of each species in the products among each reaction, SIZE: nspec*nreac
    ! integer(kind=4) :: n_ElemReac, n_ThreeBody, n_Falloff, n_PLOG    !< numbers of each type of reactions
    ! integer(kind=4) :: n_SimpleFalloff, n_TROEFalloff, n_SRI3Falloff, n_SRI5Falloff  !< numbers of each type of reactions
    ! integer(kind=4), allocatable, dimension(:) :: index_ElemReac, index_ThreeBody, index_Falloff, index_PLOG  !< index of these types of reactions 
    ! integer(kind=4), allocatable, dimension(:) :: index_SimpleFalloff, index_TROEFalloff, index_SRI3Falloff, index_SRI5Falloff    !< index of these types of reactions 
    character(len=50), parameter :: Total_reac_type = 'reac_type.dat'  !< numbers of different reaction types
    character(len=50), parameter :: reac_file_suffix = '_reac.dat'  !< reaction data file suffix of specific reaction
    real(kind=8), allocatable, dimension(:,:) :: rate_coeff             !< SIZE: nreac*11, rate parameterization coeffs for Elementary reactions, [A, b, E] in [cm, s, kmol, J/kmol], 1*3
                                                                        !< or rate parameterization coeffs for Three-body reactions, [A, b, E] in [cm, s, kmol, J/kmol], 1*3
                                                                        !< or rate parameterization coeffs for Simple Falloff reactions, [A, b, E] at low and high pressures, 1*6
                                                                        !< or rate parameterization coeffs for TROE Falloff reactions, [A, b, E] at low and high pressures and four parameters, 1*10
                                                                        !< or rate parameterization coeffs for SRI Falloff reactions, [A, b, E] at low and high pressures and three parameters, 1*11, the last two column is 1 and 0 for SRI3
    real(kind=8), allocatable, dimension(:,:) :: collision_efficiency   !< collision efficiency for Three-body reactions, defaults 1.0 for these reactions, including specific numbers from the mech, SIZE: nreac*nspec 
    integer(kind=4), allocatable, dimension(:) :: Pnum_PLOG             !< numbers of pressure points for all reactions, SIZE: nreac*1
    integer(kind=4) :: max_Pnum
    real(kind=8), allocatable, dimension(:,:,:) :: rate_coeff_PLOG      !< rate parameterization coeffs for PLOG reactions, [A, b, E] at different pressures, SIZE: nreac*max_Pnum*4

    !=========================
    ! variables for solver
    !=========================
    integer(kind=4) :: dimx_total, dimy_total, dimz_total  !< Number of grid points in each direction
    real(kind=8), allocatable, dimension(:,:,:,:) :: lisol_mass_frac    !< Mass fraction [-], nx*ny*nz*nsec
    real(kind=8), allocatable, dimension(:,:,:,:) :: lisol_mole_frac    !< Mole fraction [-], nx*ny*nz*nsec
    real(kind=8), allocatable, dimension(:,:,:) :: lisol_pressure       !< Pressure [Pa]
    real(kind=8), allocatable, dimension(:,:,:) :: lisol_temperature    !< Temperature [K]
    real(kind=8), allocatable, dimension(:) :: lisol_molmas             !< Molecular weight [kg/kmol]
    real(kind=8), allocatable, dimension(:,:,:) :: lisol_meanmolmas     !< Mean molecular weight [kg/kmol]
    real(kind=8), allocatable, dimension(:,:,:) :: lisol_density        !< Density [kg/m3]
    real(kind=8), allocatable, dimension(:,:,:,:) :: lisol_cp_species   !< Cp of each species [J/kg/K]
    real(kind=8), allocatable, dimension(:,:,:) :: lisol_cp             !< Cp of mixtures [J/kg/K]
    real(kind=8), allocatable, dimension(:,:,:) :: lisol_cv             !< Cv of mixtures [J/kg/K]
    real(kind=8), allocatable, dimension(:,:,:) :: lisol_gamma          !< specific heat [-]
    real(kind=8), allocatable, dimension(:,:,:) :: lisol_enthalpy       !< Enthalpy of mixtures [J/kg]
    real(kind=8), allocatable, dimension(:,:,:,:) :: lisol_enthalpy_spec!< Enthalpy of each species [J/kg] 
    real(kind=8), allocatable, dimension(:,:,:) :: lisol_lambda         !< Thermal conductivity [W/m]
    real(kind=8), allocatable, dimension(:,:,:) :: lisol_viscosity      !< Dynamic viscosity [Pa*s]
    real(kind=8), allocatable, dimension(:,:,:,:) :: lisol_mixDiff      !< Mixture-averaged diffusion coefficients [m^2/s]
    real(kind=8), allocatable, dimension(:,:,:) :: lisol_heat_prod_rate !< temperature production rate [W/m3]
    real(kind=8), allocatable, dimension(:,:,:) :: lisol_heat_release   !< heat release rate [W/m3]
    real(kind=8), allocatable, dimension(:,:,:,:) :: lisol_mass_prod_rate !< species mass production rate [kg/m3/s]
    real(kind=8) :: init_pressure, init_temperature
    real(kind=8), allocatable, dimension(:) :: init_frac
    character(len=500) :: init_comp, mass_or_mole, inert_gas
    logical :: ideal_gas
    integer(kind=4) :: y_forcing_type
    
    !=========================
    ! constants
    !=========================
    character(len=8), parameter :: lisol_input_file = 'LISOL_IN'
    integer(kind=4), parameter :: trans_atom      = 0
    integer(kind=4), parameter :: trans_linear    = 1
    integer(kind=4), parameter :: trans_nonlinear = 2
    real(kind=8), parameter :: SmallNumber = 1.d-300
    real(kind=8), parameter :: m2A = 1.d10  !< meter to Angstrom
    real(kind=8), parameter :: pi = 3.14159265358979323846d0 !dacos(-1.d0)  !< pi
    real(kind=8), parameter :: LightSpeed = 299792458.d0  !< light speed, [m/s]
    real(kind=8), parameter :: Avogadro = 6.02214129d26   !< Avogadro constant, [1/kmol]
    real(kind=8), parameter :: GasConstant = 8314.4621d0  !< gas constant, [J/kmol/K]
    real(kind=8), parameter :: Boltzmann = GasConstant / Avogadro  !< Boltzmann constant, [J/K]
    real(kind=8), parameter :: Debye2CoulombM = 1.d-21 / LightSpeed !< convert Debye to Coulomb*m
    real(kind=8), parameter :: OneAtm = 101325.d0 !< Pa
    ! Permeability of free space \f$ \mu_0 \f$ in N/A^2.
    real(kind=8), parameter :: permeability_0 = 4.0d-7 * pi
    ! Permittivity of free space \f$ \epsilon_0 \f$ in F/m.
    real(kind=8), parameter :: epsilon_0 = 1.0d0 / (LightSpeed * LightSpeed * permeability_0)
    real(kind=8), parameter :: cal2J = 4.184d0  !< convert calorie to J
    character(len=4), parameter :: mass_description = 'mass'  !< initial fraction based on mass
    character(len=4), parameter :: mole_description = 'mole'  !< initial fraction based on mole
    integer(kind=4), parameter :: y_forcing_sum   = 0         !< force the sum of mass fraction to unity by normalization
    integer(kind=4), parameter :: y_forcing_inert = 1         !< force the sum of mass fraction to unity through inert gas
    real(kind=8) :: omega22_table(37,8)            !< Omega22* table from Monchick and Mason, 1961
    real(kind=8) :: delta_star(8)                  !< normalized delta from Monchick and Mason, 1961
    real(kind=8) :: tstar22(37)                    !< normalized temperature for Omega22 from Monchick and Mason, 1961
    real(kind=8) :: tstar(39)                      !< extended normalized temperature for A* from Monchick and Mason, 1961
    real(kind=8) :: astar_table(39,8)              !< A* table from Monchick and Mason, 1961
    real(kind=8) :: bstar_table(39,8)              !< B* table from Monchick and Mason, 1961
    real(kind=8) :: cstar_table(39,8)              !< C* table from Monchick and Mason, 1961
    !---- the reaction types that Cantera2.4.0 supports and can be easily realized here
    integer(kind=4) :: ELEMENTARY_RXN = 100
    integer(kind=4) :: THREE_BODY_RXN = 200
    integer(kind=4) :: FALLOFF_RXN = 400
    !--three types of falloff reactions: Lindemann, TROE, and SRI
    integer(kind=4) :: SIMPLE_FALLOFF = 410
    integer(kind=4) :: TROE_FALLOFF = 420
    integer(kind=4) :: SRI_FALLOFF_3 = 433
    integer(kind=4) :: SRI_FALLOFF_5 = 435
    !--
    integer(kind=4) :: PLOG_RXN = 500
    integer(kind=4) :: non = 0
    !----
    
    !=========================
    ! others
    !=========================
    character(len=300) :: dir_comm
    character(len=300) :: dir_comm1
    ! ----- used for collision integral and subsequent transport property calculation
    integer(kind=4), parameter :: fitDelta_order = 6
    real(kind=8) :: m_logTemp(37), m_o22poly(37,fitDelta_order+1)
    real(kind=8) :: m_apoly(39,fitDelta_order+1), m_bpoly(39,fitDelta_order+1), m_cpoly(39,fitDelta_order+1)
    real(kind=8), allocatable, dimension(:,:) :: m_epsilon, m_delta, m_reducedMass, m_diam, m_dipole
    real(kind=8), allocatable, dimension(:) :: m_crot
    logical, allocatable, dimension(:) :: m_polar
    integer(kind=4), parameter :: fitTrans_order = 4
    real(kind=8), allocatable, dimension(:,:) :: m_visccoeffs, m_condcoeffs, m_diffcoeffs
    ! -----
    logical, allocatable, dimension(:,:,:) :: flag
    integer(kind=4), allocatable, dimension(:) :: devia_coeffs
    real(kind=8), allocatable, dimension(:,:,:) :: lisol_work, lisol_work2, lisol_work3, lisol_work4, lisol_work5
    real(kind=8), allocatable, dimension(:,:,:,:) :: lisol_props_species, lisol_props_species2
    real(kind=8), allocatable, dimension(:,:,:,:,:) :: m_bdiff  !< binary diffusion coefficients at unit pressure, SIZE: sub_nx*sub_ny*sub_nz*nspec*nspec

    !---- for debug
    logical :: if_reac_debug
    integer(kind=4) :: reac_debug

namelist /nml_chemistry/    mech_name, &
                            init_pressure, &
                            init_temperature, &
                            mass_or_mole, &
                            init_comp, &
                            ideal_gas, &
                            y_forcing_type, &
                            inert_gas, &
                            if_reac_debug, &
                            reac_debug

namelist /nml_initialization/   dimx_total, &
                                dimy_total, &
                                dimz_total

contains

!=============================================
! This subroutine is used for initialize the 
! parameters for collision integral calculation
!=============================================
subroutine collision_integral_para_init()
    implicit none
    real(kind=8) :: omega22_table_temp(8,37)
    real(kind=8) :: astar_table_temp(8,39), bstar_table_temp(8,39), cstar_table_temp(8,39)
    integer(kind=4) :: i, j
    real(kind=8) :: d, f_eps, f_sigma

    write(*,'(a)') '  Initial the parameters for collision integral calculation'

    omega22_table_temp  =  reshape([4.10050d0, 4.26600d0, 4.83300d0, 5.74200d0, 6.72900d0, 8.62400d0, 10.34000d0, 11.89000d0, &
                                    3.26260d0, 3.30500d0, 3.51600d0, 3.91400d0, 4.43300d0, 5.57000d0, 6.63700d0, 7.61800d0, &
                                    2.83990d0, 2.83600d0, 2.93600d0, 3.16800d0, 3.51100d0, 4.32900d0, 5.12600d0, 5.87400d0, &
                                    2.53100d0, 2.52200d0, 2.58600d0, 2.74900d0, 3.00400d0, 3.64000d0, 4.28200d0, 4.89500d0, &
                                    2.28370d0, 2.27700d0, 2.32900d0, 2.46000d0, 2.66500d0, 3.18700d0, 3.72700d0, 4.24900d0, &
                                    2.08380d0, 2.08100d0, 2.13000d0, 2.24300d0, 2.41700d0, 2.86200d0, 3.32900d0, 3.78600d0, &
                                    1.92200d0, 1.92400d0, 1.97000d0, 2.07200d0, 2.22500d0, 2.61400d0, 3.02800d0, 3.43500d0, &
                                    1.79020d0, 1.79500d0, 1.84000d0, 1.93400d0, 2.07000d0, 2.41700d0, 2.78800d0, 3.15600d0, &
                                    1.68230d0, 1.68900d0, 1.73300d0, 1.82000d0, 1.94400d0, 2.25800d0, 2.59600d0, 2.93300d0, &
                                    1.59290d0, 1.60100d0, 1.64400d0, 1.72500d0, 1.83800d0, 2.12400d0, 2.43500d0, 2.74600d0, &
                                    1.45510d0, 1.46500d0, 1.50400d0, 1.57400d0, 1.67000d0, 1.91300d0, 2.18100d0, 2.45100d0, &
                                    1.35510d0, 1.36500d0, 1.40000d0, 1.46100d0, 1.54400d0, 1.75400d0, 1.98900d0, 2.22800d0, &
                                    1.28000d0, 1.28900d0, 1.32100d0, 1.37400d0, 1.44700d0, 1.63000d0, 1.83800d0, 2.05300d0, &
                                    1.22190d0, 1.23100d0, 1.25900d0, 1.30600d0, 1.37000d0, 1.53200d0, 1.71800d0, 1.91200d0, &
                                    1.17570d0, 1.18400d0, 1.20900d0, 1.25100d0, 1.30700d0, 1.45100d0, 1.61800d0, 1.79500d0, &
                                    1.09330d0, 1.10000d0, 1.11900d0, 1.15000d0, 1.19300d0, 1.30400d0, 1.43500d0, 1.57800d0, &
                                    1.03880d0, 1.04400d0, 1.05900d0, 1.08300d0, 1.11700d0, 1.20400d0, 1.31000d0, 1.42800d0, &
                                    0.99963d0, 1.00400d0, 1.01600d0, 1.03500d0, 1.06200d0, 1.13300d0, 1.22000d0, 1.31900d0, &
                                    0.96988d0, 0.97320d0, 0.98300d0, 0.99910d0, 1.02100d0, 1.07900d0, 1.15300d0, 1.23600d0, &
                                    0.92676d0, 0.92910d0, 0.93600d0, 0.94730d0, 0.96280d0, 1.00500d0, 1.05800d0, 1.12100d0, &
                                    0.89616d0, 0.89790d0, 0.90300d0, 0.91140d0, 0.92300d0, 0.95450d0, 0.99550d0, 1.04400d0, &
                                    0.87272d0, 0.87410d0, 0.87800d0, 0.88450d0, 0.89350d0, 0.91810d0, 0.95050d0, 0.98930d0, &
                                    0.85379d0, 0.85490d0, 0.85800d0, 0.86320d0, 0.87030d0, 0.89010d0, 0.91640d0, 0.94820d0, &
                                    0.83795d0, 0.83880d0, 0.84140d0, 0.84560d0, 0.85150d0, 0.86780d0, 0.88950d0, 0.91600d0, &
                                    0.82435d0, 0.82510d0, 0.82730d0, 0.83080d0, 0.83560d0, 0.84930d0, 0.86760d0, 0.89010d0, &
                                    0.80184d0, 0.80240d0, 0.80390d0, 0.80650d0, 0.81010d0, 0.82010d0, 0.83370d0, 0.85040d0, &
                                    0.78363d0, 0.78400d0, 0.78520d0, 0.78720d0, 0.78990d0, 0.79760d0, 0.80810d0, 0.82120d0, &
                                    0.76834d0, 0.76870d0, 0.76960d0, 0.77120d0, 0.77330d0, 0.77940d0, 0.78780d0, 0.79830d0, &
                                    0.75518d0, 0.75540d0, 0.75620d0, 0.75750d0, 0.75920d0, 0.76420d0, 0.77110d0, 0.77970d0, &
                                    0.74364d0, 0.74380d0, 0.74450d0, 0.74550d0, 0.74700d0, 0.75120d0, 0.75690d0, 0.76420d0, &
                                    0.71982d0, 0.72000d0, 0.72040d0, 0.72110d0, 0.72210d0, 0.72500d0, 0.72890d0, 0.73390d0, &
                                    0.70097d0, 0.70110d0, 0.70140d0, 0.70190d0, 0.70260d0, 0.70470d0, 0.70760d0, 0.71120d0, &
                                    0.68545d0, 0.68550d0, 0.68580d0, 0.68610d0, 0.68670d0, 0.68830d0, 0.69050d0, 0.69320d0, &
                                    0.67232d0, 0.67240d0, 0.67260d0, 0.67280d0, 0.67330d0, 0.67430d0, 0.67620d0, 0.67840d0, &
                                    0.65099d0, 0.65100d0, 0.65120d0, 0.65130d0, 0.65160d0, 0.65240d0, 0.65340d0, 0.65460d0, &
                                    0.61397d0, 0.61410d0, 0.61430d0, 0.61450d0, 0.61470d0, 0.61480d0, 0.61480d0, 0.61470d0, &
                                    0.58870d0, 0.58890d0, 0.58940d0, 0.59000d0, 0.59030d0, 0.59010d0, 0.58950d0, 0.58850d0],&
                                    shape(omega22_table_temp))
    omega22_table(:,:) = transpose(omega22_table_temp(:,:))

    delta_star = reshape([0.0d0, 0.25d0, 0.50d0, 0.75d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0], shape(delta_star))

    tstar22 = reshape([0.1d0,  0.2d0,  0.3d0,  0.4d0,  0.5d0,  0.6d0,  0.7d0,  0.8d0,   0.9d0, 1.0d0, &
                       1.2d0,  1.4d0,  1.6d0,  1.8d0,  2.0d0,  2.5d0,  3.0d0,  3.5d0,   4.0d0, &
                       5.0d0,  6.0d0,  7.0d0,  8.0d0,  9.0d0, 10.0d0, 12.0d0, 14.0d0,  16.0d0, &
                      18.0d0, 20.0d0, 25.0d0, 30.0d0, 35.0d0, 40.0d0, 50.0d0, 75.0d0, 100.0d0], shape(tstar22))

    tstar = reshape([0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, 0.9d0, 1.0d0, &
                     1.2d0, 1.4d0, 1.6d0, 1.8d0, 2.0d0, 2.5d0, 3.0d0, 3.5d0, 4.0d0, &
                     5.0d0, 6.0d0, 7.0d0, 8.0d0, 9.0d0, 10.0d0, 12.0d0, 14.0d0, 16.0d0, &
                     18.0d0, 20.0d0, 25.0d0, 30.0d0, 35.0d0, 40.0d0, 50.0d0, 75.0d0, 100.0d0, 500.0d0], shape(tstar))

    astar_table_temp = reshape([1.00650d0, 1.08400d0, 1.08400d0, 1.08400d0, 1.08400d0, 1.08400d0, 1.08400d0, 1.08400d0, &
                                1.02310d0, 1.06600d0, 1.03800d0, 1.04000d0, 1.04300d0, 1.05000d0, 1.05200d0, 1.05100d0, &
                                1.04240d0, 1.04500d0, 1.04800d0, 1.05200d0, 1.05600d0, 1.06500d0, 1.06600d0, 1.06400d0, &
                                1.07190d0, 1.06700d0, 1.06000d0, 1.05500d0, 1.05800d0, 1.06800d0, 1.07100d0, 1.07100d0, &
                                1.09360d0, 1.08700d0, 1.07700d0, 1.06900d0, 1.06800d0, 1.07500d0, 1.07800d0, 1.07800d0, &
                                1.10530d0, 1.09800d0, 1.08800d0, 1.08000d0, 1.07800d0, 1.08200d0, 1.08400d0, 1.08400d0, &
                                1.11040d0, 1.10400d0, 1.09600d0, 1.08900d0, 1.08600d0, 1.08900d0, 1.09000d0, 1.09000d0, &
                                1.11140d0, 1.10700d0, 1.10000d0, 1.09500d0, 1.09300d0, 1.09500d0, 1.09600d0, 1.09500d0, &
                                1.11040d0, 1.10700d0, 1.10200d0, 1.09900d0, 1.09800d0, 1.10000d0, 1.10000d0, 1.09900d0, &
                                1.10860d0, 1.10600d0, 1.10200d0, 1.10100d0, 1.10100d0, 1.10500d0, 1.10500d0, 1.10400d0, &
                                1.10630d0, 1.10400d0, 1.10300d0, 1.10300d0, 1.10400d0, 1.10800d0, 1.10900d0, 1.10800d0, &
                                1.10200d0, 1.10200d0, 1.10300d0, 1.10500d0, 1.10700d0, 1.11200d0, 1.11500d0, 1.11500d0, &
                                1.09850d0, 1.09900d0, 1.10100d0, 1.10400d0, 1.10800d0, 1.11500d0, 1.11900d0, 1.12000d0, &
                                1.09600d0, 1.09600d0, 1.09900d0, 1.10300d0, 1.10800d0, 1.11600d0, 1.12100d0, 1.12400d0, &
                                1.09430d0, 1.09500d0, 1.09900d0, 1.10200d0, 1.10800d0, 1.11700d0, 1.12300d0, 1.12600d0, &
                                1.09340d0, 1.09400d0, 1.09700d0, 1.10200d0, 1.10700d0, 1.11600d0, 1.12300d0, 1.12800d0, &
                                1.09260d0, 1.09400d0, 1.09700d0, 1.09900d0, 1.10500d0, 1.11500d0, 1.12300d0, 1.13000d0, &
                                1.09340d0, 1.09500d0, 1.09700d0, 1.09900d0, 1.10400d0, 1.11300d0, 1.12200d0, 1.12900d0, &
                                1.09480d0, 1.09600d0, 1.09800d0, 1.10000d0, 1.10300d0, 1.11200d0, 1.11900d0, 1.12700d0, &
                                1.09650d0, 1.09700d0, 1.09900d0, 1.10100d0, 1.10400d0, 1.11000d0, 1.11800d0, 1.12600d0, &
                                1.09970d0, 1.10000d0, 1.10100d0, 1.10200d0, 1.10500d0, 1.11000d0, 1.11600d0, 1.12300d0, &
                                1.10250d0, 1.10300d0, 1.10400d0, 1.10500d0, 1.10600d0, 1.11000d0, 1.11500d0, 1.12100d0, &
                                1.10500d0, 1.10500d0, 1.10600d0, 1.10700d0, 1.10800d0, 1.11100d0, 1.11500d0, 1.12000d0, &
                                1.10720d0, 1.10700d0, 1.10800d0, 1.10800d0, 1.10900d0, 1.11200d0, 1.11500d0, 1.11900d0, &
                                1.10910d0, 1.10900d0, 1.10900d0, 1.11000d0, 1.11100d0, 1.11300d0, 1.11500d0, 1.11900d0, &
                                1.11070d0, 1.11100d0, 1.11100d0, 1.11100d0, 1.11200d0, 1.11400d0, 1.11600d0, 1.11900d0, &
                                1.11330d0, 1.11400d0, 1.11300d0, 1.11400d0, 1.11400d0, 1.11500d0, 1.11700d0, 1.11900d0, &
                                1.11540d0, 1.11500d0, 1.11600d0, 1.11600d0, 1.11600d0, 1.11700d0, 1.11800d0, 1.12000d0, &
                                1.11720d0, 1.11700d0, 1.11700d0, 1.11800d0, 1.11800d0, 1.11800d0, 1.11900d0, 1.12000d0, &
                                1.11860d0, 1.11900d0, 1.11900d0, 1.11900d0, 1.11900d0, 1.11900d0, 1.12000d0, 1.12100d0, &
                                1.11990d0, 1.12000d0, 1.12000d0, 1.12000d0, 1.12000d0, 1.12100d0, 1.12100d0, 1.12200d0, &
                                1.12230d0, 1.12200d0, 1.12200d0, 1.12200d0, 1.12200d0, 1.12300d0, 1.12300d0, 1.12400d0, &
                                1.12430d0, 1.12400d0, 1.12400d0, 1.12400d0, 1.12400d0, 1.12400d0, 1.12500d0, 1.12500d0, &
                                1.12590d0, 1.12600d0, 1.12600d0, 1.12600d0, 1.12600d0, 1.12600d0, 1.12600d0, 1.12600d0, &
                                1.12730d0, 1.12700d0, 1.12700d0, 1.12700d0, 1.12700d0, 1.12700d0, 1.12700d0, 1.12800d0, &
                                1.12970d0, 1.13000d0, 1.13000d0, 1.13000d0, 1.13000d0, 1.13000d0, 1.13000d0, 1.12900d0, &
                                1.13390d0, 1.13400d0, 1.13400d0, 1.13500d0, 1.13500d0, 1.13400d0, 1.13400d0, 1.13200d0, &
                                1.13640d0, 1.13700d0, 1.13700d0, 1.13800d0, 1.13900d0, 1.13800d0, 1.13700d0, 1.13500d0, &
                                1.14187d0, 1.14187d0, 1.14187d0, 1.14187d0, 1.14187d0, 1.14187d0, 1.14187d0, 1.14187d0],&
                                shape(astar_table_temp))
    astar_table(:,:) = transpose(astar_table_temp(:,:))
    
    bstar_table_temp = reshape([1.18520d0, 1.29630d0, 1.29630d0, 1.29630d0, 1.29630d0, 1.29630d0, 1.29630d0, 1.29630d0, &
                                1.19600d0, 1.21600d0, 1.23700d0, 1.26900d0, 1.28500d0, 1.29000d0, 1.29700d0, 1.29400d0, &
                                1.24510d0, 1.25700d0, 1.34000d0, 1.38900d0, 1.36600d0, 1.32700d0, 1.31400d0, 1.27800d0, &
                                1.29000d0, 1.29400d0, 1.27200d0, 1.25800d0, 1.26200d0, 1.28200d0, 1.29000d0, 1.29900d0, &
                                1.29860d0, 1.29100d0, 1.28400d0, 1.27800d0, 1.27700d0, 1.28800d0, 1.29400d0, 1.29700d0, &
                                1.28650d0, 1.28100d0, 1.27600d0, 1.27200d0, 1.27700d0, 1.28600d0, 1.29200d0, 1.29800d0, &
                                1.26650d0, 1.26400d0, 1.26100d0, 1.26300d0, 1.26900d0, 1.28400d0, 1.29200d0, 1.29800d0, &
                                1.24550d0, 1.24400d0, 1.24800d0, 1.25500d0, 1.26200d0, 1.27800d0, 1.28900d0, 1.29600d0, &
                                1.22530d0, 1.22500d0, 1.23400d0, 1.24000d0, 1.25200d0, 1.27100d0, 1.28400d0, 1.29500d0, &
                                1.20780d0, 1.21000d0, 1.21600d0, 1.22700d0, 1.24200d0, 1.26400d0, 1.28100d0, 1.29200d0, &
                                1.19190d0, 1.19200d0, 1.20500d0, 1.21600d0, 1.23000d0, 1.25600d0, 1.27300d0, 1.28700d0, &
                                1.16780d0, 1.17200d0, 1.18100d0, 1.19500d0, 1.20900d0, 1.23700d0, 1.26100d0, 1.27700d0, &
                                1.14960d0, 1.15500d0, 1.16100d0, 1.17400d0, 1.18900d0, 1.22100d0, 1.24600d0, 1.26600d0, &
                                1.13660d0, 1.14100d0, 1.14700d0, 1.15900d0, 1.17400d0, 1.20200d0, 1.23100d0, 1.25600d0, &
                                1.12700d0, 1.13000d0, 1.13800d0, 1.14800d0, 1.16200d0, 1.19100d0, 1.21800d0, 1.24200d0, &
                                1.11970d0, 1.12200d0, 1.12900d0, 1.14000d0, 1.14900d0, 1.17800d0, 1.20500d0, 1.23100d0, &
                                1.10800d0, 1.11000d0, 1.11600d0, 1.12200d0, 1.13200d0, 1.15400d0, 1.18000d0, 1.20500d0, &
                                1.10160d0, 1.10300d0, 1.10700d0, 1.11200d0, 1.12000d0, 1.13800d0, 1.16000d0, 1.18300d0, &
                                1.09800d0, 1.09900d0, 1.10200d0, 1.10600d0, 1.11200d0, 1.12700d0, 1.14500d0, 1.16500d0, &
                                1.09580d0, 1.09700d0, 1.09900d0, 1.10200d0, 1.10700d0, 1.11900d0, 1.13500d0, 1.15300d0, &
                                1.09350d0, 1.09400d0, 1.09500d0, 1.09700d0, 1.10000d0, 1.10900d0, 1.12000d0, 1.13400d0, &
                                1.09250d0, 1.09200d0, 1.09400d0, 1.09500d0, 1.09800d0, 1.10400d0, 1.11200d0, 1.12200d0, &
                                1.09220d0, 1.09200d0, 1.09300d0, 1.09400d0, 1.09600d0, 1.10000d0, 1.10600d0, 1.11500d0, &
                                1.09220d0, 1.09200d0, 1.09300d0, 1.09300d0, 1.09500d0, 1.09800d0, 1.10300d0, 1.11000d0, &
                                1.09230d0, 1.09200d0, 1.09300d0, 1.09300d0, 1.09400d0, 1.09700d0, 1.10100d0, 1.10600d0, &
                                1.09230d0, 1.09200d0, 1.09200d0, 1.09300d0, 1.09400d0, 1.09600d0, 1.09900d0, 1.10300d0, &
                                1.09270d0, 1.09300d0, 1.09300d0, 1.09300d0, 1.09400d0, 1.09500d0, 1.09800d0, 1.10100d0, &
                                1.09300d0, 1.09300d0, 1.09300d0, 1.09300d0, 1.09400d0, 1.09400d0, 1.09600d0, 1.09900d0, &
                                1.09330d0, 1.09400d0, 1.09300d0, 1.09400d0, 1.09400d0, 1.09500d0, 1.09600d0, 1.09800d0, &
                                1.09370d0, 1.09300d0, 1.09400d0, 1.09400d0, 1.09400d0, 1.09400d0, 1.09600d0, 1.09700d0, &
                                1.09390d0, 1.09400d0, 1.09400d0, 1.09400d0, 1.09400d0, 1.09500d0, 1.09500d0, 1.09700d0, &
                                1.09430d0, 1.09400d0, 1.09400d0, 1.09400d0, 1.09500d0, 1.09500d0, 1.09600d0, 1.09600d0, &
                                1.09440d0, 1.09500d0, 1.09400d0, 1.09400d0, 1.09400d0, 1.09500d0, 1.09500d0, 1.09600d0, &
                                1.09440d0, 1.09400d0, 1.09500d0, 1.09400d0, 1.09400d0, 1.09500d0, 1.09600d0, 1.09600d0, &
                                1.09430d0, 1.09500d0, 1.09400d0, 1.09400d0, 1.09500d0, 1.09500d0, 1.09500d0, 1.09500d0, &
                                1.09410d0, 1.09400d0, 1.09400d0, 1.09400d0, 1.09400d0, 1.09400d0, 1.09400d0, 1.09600d0, &
                                1.09470d0, 1.09500d0, 1.09400d0, 1.09400d0, 1.09300d0, 1.09300d0, 1.09400d0, 1.09500d0, &
                                1.09570d0, 1.09500d0, 1.09400d0, 1.09300d0, 1.09200d0, 1.09300d0, 1.09300d0, 1.09400d0, &
                                1.10185d0, 1.10185d0, 1.10185d0, 1.10185d0, 1.10185d0, 1.10185d0, 1.10185d0, 1.10185d0],&
                                shape(bstar_table_temp))
    bstar_table(:,:) = transpose(bstar_table_temp(:,:))
    
    cstar_table_temp = reshape([0.88890d0, 0.77778d0, 0.77778d0, 0.77778d0, 0.77778d0, 0.77778d0, 0.77778d0, 0.77778d0, &
                                0.88575d0, 0.89880d0, 0.83780d0, 0.80290d0, 0.78760d0, 0.78050d0, 0.77990d0, 0.78010d0, &
                                0.87268d0, 0.86920d0, 0.86470d0, 0.84790d0, 0.82370d0, 0.79750d0, 0.78810d0, 0.77840d0, &
                                0.85182d0, 0.85250d0, 0.83660d0, 0.81980d0, 0.80540d0, 0.79030d0, 0.78390d0, 0.78200d0, &
                                0.83542d0, 0.83620d0, 0.83060d0, 0.81960d0, 0.80760d0, 0.79180d0, 0.78420d0, 0.78060d0, &
                                0.82629d0, 0.82780d0, 0.82520d0, 0.81690d0, 0.80740d0, 0.79160d0, 0.78380d0, 0.78020d0, &
                                0.82299d0, 0.82490d0, 0.82300d0, 0.81650d0, 0.80720d0, 0.79220d0, 0.78390d0, 0.77980d0, &
                                0.82357d0, 0.82570d0, 0.82410d0, 0.81780d0, 0.80840d0, 0.79270d0, 0.78390d0, 0.77940d0, &
                                0.82657d0, 0.82800d0, 0.82640d0, 0.81990d0, 0.81070d0, 0.79390d0, 0.78420d0, 0.77960d0, &
                                0.83110d0, 0.82340d0, 0.82950d0, 0.82280d0, 0.81360d0, 0.79600d0, 0.78540d0, 0.77980d0, &
                                0.83630d0, 0.83660d0, 0.83420d0, 0.82670d0, 0.81680d0, 0.79860d0, 0.78640d0, 0.78050d0, &
                                0.84762d0, 0.84740d0, 0.84380d0, 0.83580d0, 0.82500d0, 0.80410d0, 0.79040d0, 0.78220d0, &
                                0.85846d0, 0.85830d0, 0.85300d0, 0.84440d0, 0.83360d0, 0.81180d0, 0.79570d0, 0.78540d0, &
                                0.86840d0, 0.86740d0, 0.86190d0, 0.85310d0, 0.84230d0, 0.81860d0, 0.80110d0, 0.78980d0, &
                                0.87713d0, 0.87550d0, 0.87090d0, 0.86160d0, 0.85040d0, 0.82650d0, 0.80720d0, 0.79390d0, &
                                0.88479d0, 0.88310d0, 0.87790d0, 0.86950d0, 0.85780d0, 0.83380d0, 0.81330d0, 0.79900d0, &
                                0.89972d0, 0.89860d0, 0.89360d0, 0.88460d0, 0.87420d0, 0.85040d0, 0.82940d0, 0.81250d0, &
                                0.91028d0, 0.90890d0, 0.90430d0, 0.89670d0, 0.88690d0, 0.86490d0, 0.84380d0, 0.82530d0, &
                                0.91793d0, 0.91660d0, 0.91250d0, 0.90580d0, 0.89700d0, 0.87680d0, 0.85570d0, 0.83720d0, &
                                0.92371d0, 0.92260d0, 0.91890d0, 0.91280d0, 0.90500d0, 0.88610d0, 0.86640d0, 0.84840d0, &
                                0.93135d0, 0.93040d0, 0.92740d0, 0.92260d0, 0.91640d0, 0.90060d0, 0.88330d0, 0.86620d0, &
                                0.93607d0, 0.93530d0, 0.93290d0, 0.92910d0, 0.92400d0, 0.91090d0, 0.89580d0, 0.88020d0, &
                                0.93927d0, 0.93870d0, 0.93660d0, 0.93340d0, 0.92920d0, 0.91620d0, 0.90500d0, 0.89110d0, &
                                0.94149d0, 0.94090d0, 0.93930d0, 0.93660d0, 0.93310d0, 0.92360d0, 0.91220d0, 0.89970d0, &
                                0.94306d0, 0.94260d0, 0.94120d0, 0.93880d0, 0.93570d0, 0.92760d0, 0.91750d0, 0.90650d0, &
                                0.94419d0, 0.94370d0, 0.94250d0, 0.94060d0, 0.93800d0, 0.93080d0, 0.92190d0, 0.91190d0, &
                                0.94571d0, 0.94550d0, 0.94450d0, 0.94300d0, 0.94090d0, 0.93530d0, 0.92830d0, 0.92010d0, &
                                0.94662d0, 0.94640d0, 0.94560d0, 0.94440d0, 0.94280d0, 0.93820d0, 0.93250d0, 0.92580d0, &
                                0.94723d0, 0.94710d0, 0.94640d0, 0.94550d0, 0.94420d0, 0.94050d0, 0.93550d0, 0.92980d0, &
                                0.94764d0, 0.94740d0, 0.94690d0, 0.94620d0, 0.94500d0, 0.94180d0, 0.93780d0, 0.93280d0, &
                                0.94790d0, 0.94780d0, 0.94740d0, 0.94650d0, 0.94570d0, 0.94300d0, 0.93940d0, 0.93520d0, &
                                0.94827d0, 0.94810d0, 0.94800d0, 0.94720d0, 0.94670d0, 0.94470d0, 0.94220d0, 0.93910d0, &
                                0.94842d0, 0.94840d0, 0.94810d0, 0.94780d0, 0.94720d0, 0.94580d0, 0.94370d0, 0.94150d0, &
                                0.94852d0, 0.94840d0, 0.94830d0, 0.94800d0, 0.94750d0, 0.94650d0, 0.94490d0, 0.94300d0, &
                                0.94861d0, 0.94870d0, 0.94840d0, 0.94810d0, 0.94790d0, 0.94680d0, 0.94550d0, 0.94300d0, &
                                0.94872d0, 0.94860d0, 0.94860d0, 0.94830d0, 0.94820d0, 0.94750d0, 0.94640d0, 0.94520d0, &
                                0.94881d0, 0.94880d0, 0.94890d0, 0.94900d0, 0.94870d0, 0.94820d0, 0.94760d0, 0.94680d0, &
                                0.94863d0, 0.94870d0, 0.94890d0, 0.94910d0, 0.94930d0, 0.94910d0, 0.94830d0, 0.94760d0, &
                                0.94444d0, 0.94444d0, 0.94444d0, 0.94444d0, 0.94444d0, 0.94444d0, 0.94444d0, 0.94444d0],&
                                shape(cstar_table_temp))
    cstar_table(:,:) = transpose(cstar_table_temp(:,:))

    do i = 1,nspec
        m_dipole(i,i) = spec_trans_ISO(i,4)
        m_polar(i) = (spec_trans_ISO(i,4) .gt. 0.d0)
        if(spec_trans(i,1) == 0)then ! atom
            m_crot(i) = 0.d0
        elseif(spec_trans(i,1) == 1)then ! linear
            m_crot(i) = 1.d0
        elseif(spec_trans(i,1) == 2)then ! nonlinear
            m_crot(i) = 1.5d0
        endif
    enddo
    do i = 1,nspec
        do j = 1,nspec
            ! the reduced mass
            m_reducedMass(i,j) = lisol_molmas(i) * lisol_molmas(j) &
                                 / (Avogadro * (lisol_molmas(i) + lisol_molmas(j)))
            ! hard-sphere diameter for (i,j) collisions
            m_diam(i,j) = 0.5d0 * (spec_trans_ISO(i,3) + spec_trans_ISO(j,3))

            ! the effective well depth for (i,j) collisions
            m_epsilon(i,j) = dsqrt(spec_trans_ISO(i,2) * spec_trans_ISO(j,2))

            ! the effective dipole moment for (i,j) collisions
            m_dipole(i,j) = dsqrt(m_dipole(i,i) * m_dipole(j,j))

            ! reduced dipole moment delta* (nondimensional)
            d = m_diam(i,j)
            m_delta(i,j) = 0.5d0 * m_dipole(i,j) * m_dipole(i,j) &
                           / (4.d0 * pi * epsilon_0 * m_epsilon(i,j) * (d**3.d0))
            call makePolarCorrections(i, j, f_eps, f_sigma)
            m_diam(i,j) = m_diam(i,j) * f_sigma
            m_epsilon(i,j) = m_epsilon(i,j) * f_eps

            ! properties are symmetric
            m_reducedMass(j,i) = m_reducedMass(i,j)
            m_diam(j,i) = m_diam(i,j)
            m_epsilon(j,i) = m_epsilon(i,j)
            m_dipole(j,i) = m_dipole(i,j)
            m_delta(j,i) = m_delta(i,j)
        enddo
    enddo


    return
end subroutine collision_integral_para_init
!=============================================
! This subroutine is used for calculate the
! correction for cases that one molecule is polar
! and another one is non-polar
!=============================================
subroutine makePolarCorrections(i, j, f_eps, f_sigma)
    implicit none
    integer(kind=4), intent(in) :: i, j
    real(kind=8), intent(out) :: f_eps, f_sigma
    integer(kind=4) :: kp, knp
    real(kind=8) :: d3np, d3p, alpha_star, mu_p_star, xi

    ! write(*,'(a)') '  Make polar correction for well depth and diameter'

    ! no correction if both are nonpolar, or both are polar
    if(m_polar(i) .eqv. m_polar(j))then
        f_eps = 1.d0
        f_sigma = 1.d0
        return
    endif
    ! corrections to the effective diameter and well depth
    ! if one is polar and one is non-polar
    if(m_polar(i))then ! the polar one
        kp = i
    else
        kp = j
    endif
    if(i .eq. kp)then ! the nonpolar one
        knp = j
    else
        knp = i
    endif
    d3np = spec_trans_ISO(knp, 3)**3.d0
    d3p = spec_trans_ISO(kp, 3)**3.d0
    alpha_star = spec_trans_ISO(knp, 5) / d3np
    mu_p_star = m_dipole(kp,kp) / dsqrt(4.d0 * pi * epsilon_0 * d3p * spec_trans_ISO(kp, 2))
    xi = 1.d0 + 0.25d0 * alpha_star * mu_p_star * mu_p_star &
         * dsqrt(spec_trans_ISO(kp, 2) / spec_trans_ISO(knp, 2))
    f_sigma = xi**(-1.0d0/6.0d0)
    f_eps = xi * xi

    return
end subroutine makePolarCorrections
!=============================================
! This subroutine is used for setup the collision
! integral, and call the subroutine fitProperties
! to fit transport properties within the minTemp
! and maxTemp valid for thermodynamic
!=============================================
subroutine setupCollisionIntegral()
    implicit none
    ! real(kind=8) :: tstar_min = 1.0d8, tstar_max = 0.0d0, 
    real(kind=8) :: minTemp, maxTemp
    integer(kind=4) :: i, j

    write(*,'(a)') '  Setup the collision integral fit and interpolation'

    minTemp = spec_therm(1,1)
    maxTemp = spec_therm(1,3)
    do i = 2,nspec
        minTemp = dmin1(minTemp, spec_therm(i, 1))
        maxTemp = dmin1(maxTemp, spec_therm(i, 3))
    enddo
    ! write(*,*) 'maxTemp is ', maxTemp, ', minTemp is ', minTemp

    ! do i = 1,nspec
    !     do j = i,nspec
    !         tstar_min = dmin1(tstar_min, Boltzmann * minTemp / m_epsilon(i,j))
    !         tstar_max = dmax1(tstar_max, Boltzmann * maxTemp / m_epsilon(i,j))
    !     enddo
    ! enddo
    ! initialize the collision integral calculator for the desired T* range
    call collision_integral_init()
    call fitProperties(minTemp, maxTemp)

    return
end subroutine setupCollisionIntegral
!=============================================
! This subroutine is used for fitting transport 
! properties within the minTemp and maxTemp 
! valid for thermodynamic
!=============================================
subroutine fitProperties(minTemp, maxTemp)
    implicit none
    real(kind=8), intent(in) :: minTemp, maxTemp
    integer(kind=4) :: np, n, k, i
    real(kind=8) :: dt, t, visc, tstar, fz_298, sqrt_T, om22, om11, diffcoeff, cond, cp_R
    real(kind=8) :: f_int, cv_rot, A_factor, fz_tstar, B_factor, c1, cv_int, f_rot, f_trans
    real(kind=8), allocatable, dimension(:) :: tlog, spvisc, spcond, w, w2, c, c2
    logical :: flag
    real(kind=8) :: temp

    write(*,'(a)') '  Fit the transport properties'

    np = 50
    allocate(tlog(np), spvisc(np), spcond(np), w(np), w2(np))
    allocate(c(fitTrans_order+1), c2(fitTrans_order+1))
    dt = (maxTemp - minTemp) / dble(np-1)
    do n = 0, np-1
        t = minTemp + dble(dt*n)
        tlog(n+1) = dlog(t)
    enddo
    ! write(*,*) 'tlog(1) = ', tlog(1), ', tlog(np) = ', tlog(np)
    do k = 1,nspec
        tstar = Boltzmann * 298.d0 / spec_trans_ISO(k, 2)
        ! if(spec_name(k) == 'H2O') write(*,*) 'tstar_298 is ', tstar
        ! Scaling factor for temperature dependence of z_rot. [Kee2003] Eq.
        ! 12.112 or [Kee2017] Eq. 11.115
        fz_298 = 1.d0 + pi**1.5d0 / dsqrt(tstar) * (0.5d0 + 1.d0 / tstar) &
                 + (0.25d0 * pi * pi + 2.d0) / tstar
        do n = 0,np-1
            t = minTemp + dble(dt*n)
            ! get the cp_R for each species here, it's convenient in Cantera, but here we 
            ! just repeat the code in subroutine lisol_compute_cp_species for species k
            flag = (t <= spec_therm(k,2))
            cp_R =  merge((spec_therm(k, 4) + spec_therm(k, 5)*t + &
                           spec_therm(k, 6)*(t**2.d0) + spec_therm(k, 7)*(t**3.d0) + &
                           spec_therm(k, 8)*(t**4.d0)), &
                          (spec_therm(k,11) + spec_therm(k,12)*t + &
                           spec_therm(k,13)*(t**2.d0) + spec_therm(k,14)*(t**3.d0) + &
                           spec_therm(k,15)*(t**4.d0)), flag)
            ! if(spec_name(k) == 'H2O' .and. n==0) write(*,*) 'Cp_ref of H2O is ', cp_R
            tstar = Boltzmann * t / spec_trans_ISO(k, 2)
            sqrt_T = dsqrt(t)
            call omega22(tstar, m_delta(k,k), om22)
            call astar(tstar, m_delta(k,k), temp)
            om11 = om22 / temp

            ! self-diffusion coefficient, without polar corrections
            diffcoeff = 3.d0 / 16.d0 * dsqrt(2.d0 * pi / m_reducedMass(k,k)) &
                        * ((Boltzmann * t)**1.5d0) / (pi * spec_trans_ISO(k,3) * spec_trans_ISO(k,3) * om11)
            ! viscosity
            visc = 5.d0 / 16.d0 * dsqrt(pi * lisol_molmas(k) * Boltzmann * t / Avogadro) &
                   / (om22 * pi * spec_trans_ISO(k,3) * spec_trans_ISO(k,3))
            ! thermal conductivity
            f_int = lisol_molmas(k) / (GasConstant * t) * diffcoeff / visc
            cv_rot = m_crot(k)
            A_factor = 2.5d0 - f_int
            fz_tstar = 1.0d0 + pi**1.5d0 / dsqrt(tstar) * (0.5d0 + 1.d0 / tstar) &
                       + (0.25d0 * pi * pi + 2.d0) / tstar
            B_factor = spec_trans_ISO(k,6) * fz_298 / fz_tstar + 2.d0/pi * (5.d0/3.d0 * cv_rot + f_int)
            c1 = 2.d0/pi * A_factor / B_factor
            cv_int = cp_R - 2.5d0 - cv_rot
            f_rot   = f_int * (1.d0 + c1)
            f_trans = 2.5d0 * (1.d0 - c1 * cv_rot / 1.5d0)
            cond = (visc / lisol_molmas(k)) * GasConstant * (f_trans * 1.5d0 &
                                                             + f_rot * cv_rot + f_int * cv_int)
            ! the viscosity should be proportional approximately to
            ! sqrt(T); therefore, visc/sqrt(T) should have only a weak
            ! temperature dependence. And since the mixture rule requires
            ! the square root of the pure-species viscosity, fit the square
            ! root of (visc/sqrt(T)) to avoid having to compute square
            ! roots in the mixture rule.
            spvisc(n+1) = dsqrt(visc / sqrt_T)
            ! the pure-species conductivity scales approximately with
            ! sqrt(T). Unlike the viscosity, there is no reason here to fit
            ! the square root, since a different mixture rule is used.
            spcond(n+1) = cond / sqrt_T
            w(n+1) = 1.d0 / (spvisc(n+1) * spvisc(n+1))
            w2(n+1) = 1.d0 / (spcond(n+1) * spcond(n+1))
        enddo
        call lisol_polyfit(tlog(:), spvisc(:), w(:),  np, fitTrans_order, c(:))
        ! ! debug
        ! if(spec_name(k) == 'H2O')then
        !     do i = 1,fitTrans_order+1
        !         write(*,*) 'Coeffs(', i, ') for mu of H2O is ', c(i)
        !     enddo
        ! endif
        call lisol_polyfit(tlog(:), spcond(:), w2(:), np, fitTrans_order, c2(:))
        ! ! debug
        ! if(spec_name(k) == 'H2O')then
        !     do i = 1,fitTrans_order+1
        !         write(*,*) 'Coeffs(', i, ') for lambda of H2O is ', c2(i)
        !     enddo
        ! endif
        m_visccoeffs(k,:) = c(:)
        m_condcoeffs(k,:) = c2(:)
    enddo
    
    call fitDiffCoeffs(minTemp, maxTemp)

    deallocate(tlog, spvisc, spcond, w, w2, c, c2)

    return
end subroutine fitProperties
!=============================================
! This subroutine is used for fitting transport 
! properties within the minTemp and maxTemp 
! valid for thermodynamic
!=============================================
subroutine fitDiffCoeffs(minTemp, maxTemp)
    implicit none
    real(kind=8), intent(in) :: minTemp, maxTemp
    integer(kind=4) :: np, n, k, j, icount
    real(kind=8) :: dt, t, eps, tstar, sigma, om11, temp, diffcoeff
    real(kind=8), allocatable, dimension(:) :: tlog_diff, w_diff, c_diff, diff

    np = 50
    allocate(tlog_diff(np), w_diff(np))
    allocate(c_diff(fitTrans_order+1))
    allocate(diff(np))
    dt = (maxTemp - minTemp) / dble(np-1)
    do n = 0, np-1
        t = minTemp + dble(dt*n)
        tlog_diff(n+1) = dlog(t)
    enddo

    icount = 1
    do k = 1,nspec
        do j = k,nspec
            do n = 0,np-1
                t = minTemp + dble(dt*n)
                eps = m_epsilon(j,k)
                tstar = Boltzmann * t / eps
                sigma = m_diam(j,k)

                call omega22(tstar, m_delta(j,k), om11)
                call astar(tstar, m_delta(j,k), temp)
                om11 = om11 / temp

                diffcoeff = 3.d0/16.d0 * dsqrt(2.d0 * pi / m_reducedMass(k,j)) &
                            * (Boltzmann * t)**(1.5d0) / (pi * sigma * sigma * om11)

                diff(n+1) = diffcoeff / (t**(1.5d0))
                w_diff(n+1) = 1.d0 / (diff(n+1)*diff(n+1))
            enddo
            call lisol_polyfit(tlog_diff(:), diff(:), w_diff(:), np, fitTrans_order, c_diff(:))
            m_diffcoeffs(icount,:) = c_diff(:)
            icount = icount + 1
        enddo
    enddo
    
    deallocate(tlog_diff, w_diff, c_diff, diff)
    return
end subroutine fitDiffCoeffs
!=============================================
! This subroutine is used for fit the polynomial
! with given x, y, and order, used for collision
! integral fit, using column_pivot_gauss method
! NOTE: the coefficients obtained are listed from 
! lower order to higher order !!!!! (Mygo)
!=============================================
subroutine lisol_polyfit(x, y, w, n, order, p)
    ! use blas77
	! use lapack77
    implicit none
    integer(kind=4), intent(in) :: n, order
    real(kind=8), dimension(:), intent(in) :: x(n), y(n), w(n) ! weights
    real(kind=8), dimension(:), intent(out) :: p(order+1)
    integer(kind=4) :: i, j, k, max_row, ierr
    real(kind=8), allocatable, dimension(:,:) :: A, An, At
    real(kind=8), allocatable, dimension(:) :: temp_array, b, weights, ipiv
    real(kind=8) :: max_val, factor, temp

    allocate(b(order+1), weights(n), ipiv(order+1))
    allocate(A(n, order+1), At(order+1, n), temp_array(order+1))
    allocate(An(order+1, order+1))

    ! ! SQRT weights
    ! weights(:) = dsqrt(w(:))
    weights(:) = w(:)

    A(:,1) = 1.d0
    do i = 2,order+1
        A(:,i) = A(:,i-1) * x(:)
    enddo
    At(:,:) = transpose(A(:,:))
    do i = 1,n
        A(i,:) = A(i,:) * weights(i)  ! weights
    enddo
    An(:,:) = matmul(At(:,:), A(:,:))

    do i = 1,order+1
        temp = 0.d0
        do j = 1,n
            temp = temp + At(i,j) * (weights(j) * y(j))  ! weights
        enddo
        b(i) = temp
    enddo

    ! p(:) = b(:)
    ! call dgesv(order+1, 1, An, order+1, ipiv, p, order+1, ierr)
    ! if(ierr .ne. 0)then
    !     call lisol_error(ierr, 'lisol_polyfit: dgesv failed!')
    ! endif

    !> Forward elimination
    do i = 1,order
        !> 1. find the pivot
        max_row = i
        max_val = dabs(An(i, i))
        do k = i+1,order+1
            if(dabs(An(k, i)) .gt. max_val)then
                max_row = k
                max_val = dabs(An(k, i))
            endif
        enddo
        !> 2. exchange the rows
        if(max_row .ne. i)then
            temp_array(:) = An(i,:)
            An(i,:) = An(max_row,:)
            An(max_row,:) = temp_array(:)

            temp = b(i)
            b(i) = b(max_row)
            b(max_row) = temp
        endif
        !> 3. elimination
        do j = i+1,order+1
            factor = An(j, i) / An(i, i)
            An(j, i:(order+1)) = An(j, i:(order+1)) - factor * An(i, i:(order+1))
            b(j) = b(j) - factor * b(i)
        enddo
    enddo

    !> Back substitution
    p(order+1) = b(order+1) / An(order+1, order+1)
    do i = order, 1, -1
        p(i) = (b(i) - sum(An(i, (i+1):(order+1)) * p((i+1):(order+1)))) / An(i,i)
    enddo

    deallocate(b, ipiv)
    deallocate(A, At, temp_array)
    deallocate(An)  
    return
end subroutine lisol_polyfit
!=============================================
! This subroutine is used for quadratic 
! interpolation, which is used to determine the
! normalized temperature for properties fit
!=============================================
subroutine quadInterp(x0, x, y, interp)
    implicit none
    real(kind=8), intent(in) :: x0
    real(kind=8), intent(in) :: x(3), y(3)
    real(kind=8), intent(out) :: interp
    real(kind=8) :: dx21, dx32, dx31
    real(kind=8) :: dy32, dy21, a

    dx21 = x(2) - x(1)
    dx32 = x(3) - x(2)
    dx31 = dx21 + dx32
    dy32 = y(3) - y(2)
    dy21 = y(2) - y(1)
    a = (dx21*dy32 - dy21*dx32)/(dx21*dx31*dx32)
    interp = a*(x0 - x(1))*(x0 - x(2)) + (dy21/dx21)*(x0 - x(2)) + y(2)

    return
end subroutine quadInterp
!=============================================
! This subroutine is used for evaluating polynomial
! results with given order
! NOTE: the coefficients are listed from 
! lower order to higher order !!!!! (Mygo)
!=============================================
subroutine polyval(x0, c, order, y0)
    implicit none
    real(kind=8), intent(in) :: x0
    integer(kind=4), intent(in) :: order
    real(kind=8), intent(in) :: c(order+1)
    real(kind=8), intent(out) :: y0

    integer(kind=4) :: i

    y0 = 0.d0
    do i = 1,order+1
        y0 = y0 + c(i) * (x0**(dble(i-1)))
    enddo

    return
end subroutine polyval
!=============================================
! This subroutine is used for initialize the 
! collision integral calculation, including 
! 1. range of normalized temperature
! 2. fit of tables along Delta*
!=============================================
subroutine collision_integral_init()
    implicit none
    integer(kind=4) :: i
    real(kind=8) :: c(fitDelta_order+1)

    write(*,'(a)') '  Fit the integral table along normalized dipole moment'
    ! write(*,*) 'breakpoint here'
    do i = 1,37
        m_logTemp(i) = dlog(tstar(i+1))
        ! omega22
        call fitDelta(0, i, fitDelta_order, c(:))
        m_o22poly(i,:) = c(:)
    enddo
    do i = 1,39
        ! A*
        call fitDelta(1, i, fitDelta_order, c(:))
        m_apoly(i,:) = c(:)
        ! B*
        call fitDelta(2, i, fitDelta_order, c(:))
        m_bpoly(i,:) = c(:)
        ! C*
        call fitDelta(3, i, fitDelta_order, c(:))
        m_cpoly(i,:) = c(:)
    enddo

    return
end subroutine collision_integral_init
!=============================================
! This subroutine is used for fitting the table
! along Delta*, mainly determining the row index
!=============================================
subroutine fitDelta(table, ntstar, degree, c)
    implicit none
    integer(kind=4), intent(in) :: table, ntstar, degree
    real(kind=8), intent(out) :: c(degree+1)
    real(kind=8) :: w(8), y(8), x(8)

    select case(table)
        case(0)
            y(:) = omega22_table(ntstar, :)
        case(1)
            y(:) = astar_table(ntstar, :)
        case(2)
            y(:) = bstar_table(ntstar, :)
        case(3)
            y(:) = cstar_table(ntstar, :)
        case default
            call lisol_error(1, 'fitDelta: error case table!')
    end select

    w(:) = 1.d0
    x(:) = delta_star(:)
    call lisol_polyfit(delta_star(:), y(:), w(:), 8, fitDelta_order, c(:))

    return
end subroutine fitDelta
!=============================================
! This subroutine is used for getting the Omega22*
! quadratic interpolation along normalized temperature
!=============================================
subroutine omega22(ts, deltastar, omega)
    implicit none
    real(kind=8), intent(in) :: ts, deltastar
    real(kind=8), intent(out) :: omega
    integer(kind=4) :: i, i1, i2
    real(kind=8) :: values(3), temp, c(fitDelta_order+1), x(3)

    do i = 1,37
        if(ts .lt. tstar22(i)) exit
    enddo
    i1 = max(i-1, 0)
    i2 = i1 + 3
    if(i2 > 37)then
        i2 = 37
        i1 = i2 - 3
    endif
    do i = i1, i2-1
        if(deltastar .eq. 0.d0)then
            values(i-i1+1) = omega22_table(i,1)
        else
            c(:) = m_o22poly(i,:)
            call polyval(deltastar, c(:), fitDelta_order, temp)
            values(i-i1+1) = temp
        endif
    enddo

    x(:) = m_logTemp(i1:(i1+2))
    call quadInterp(dlog(ts), x(:), values(:), temp)
    omega = temp

    return
end subroutine omega22
!=============================================
! This subroutine is used for getting the A*
! quadratic interpolation along normalized temperature
!=============================================
subroutine astar(ts, deltastar, a)
    implicit none
    real(kind=8), intent(in) :: ts, deltastar
    real(kind=8), intent(out) :: a
    integer(kind=4) :: i, i1, i2
    real(kind=8) :: values(3), temp, c(fitDelta_order+1), x(3)

    do i = 1,37
        if(ts .lt. tstar22(i)) exit
    enddo
    i1 = max(i-1, 0)
    i2 = i1 + 3
    if(i2 > 37)then
        i2 = 37
        i1 = i2 - 3
    endif
    do i = i1, i2-1
        if(deltastar .eq. 0.d0)then
            values(i-i1+1) = astar_table(i+1,1)
        else
            c(:) = m_apoly(i+1,:)
            call polyval(deltastar, c(:), fitDelta_order, temp)
            values(i-i1+1) = temp
        endif
    enddo

    x(:) = m_logTemp(i1:(i1+2))
    call quadInterp(dlog(ts), x(:), values(:), temp)
    a = temp

    return
end subroutine astar
!=============================================
! This subroutine is used for print error massage and stop the program
!=============================================
subroutine lisol_error(errid, err_msg)
    implicit none
    integer(kind=4), intent(in) :: errid
    character(len=*), intent(in) :: err_msg

    write(*,*)
    write(*,*)
    write(*,'(a)') 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
    write(*,'(a)') ' --> Lisol ERROR in '//err_msg
    write(*,'(a)') 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
    write(*,*)
    write(*,'(a)') ' error code: ', errid
    write(*,*)
    write(*,'(a)') '========================================'
    write(*,'(a)') '      error termination, stopping'
    write(*,'(a)') '========================================'
    write(*,*)
    write(*,*)

    stop
    return
end subroutine lisol_error

end module lisol_chemistry_mod