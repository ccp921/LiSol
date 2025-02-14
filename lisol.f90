!=============================================
! This program test for chemistry mechanism
!=============================================
program lisol
    use cantera
    use lisol_chemistry_mod
    implicit none
    character(len=500) :: mech_file, mech_file_length, file_temp, ctemp
    integer(kind=4) :: interface_file = 921, input_file = 129
    integer(kind=4) :: i, j, l, ierr, itemp
    real(kind=8) :: temp
    
    ! read input file
    open(input_file,file=trim(lisol_input_file),status='old',position='rewind',action='read',iostat=ierr)
    if (ierr .ne. 0) then
        call lisol_error(ierr, 'lisol: Could not open input file')
    endif

    read(input_file,nml=nml_chemistry,iostat=ierr)
    rewind(input_file)
    if(ierr .ne. 0) then
        call lisol_error(ierr, 'lisol: error during read of nml_chemistry')
    endif

    read(input_file,nml=nml_initialization,iostat=ierr)
    rewind(input_file)
    if(ierr .ne. 0) then
        call lisol_error(ierr, 'lisol: error during read of nml_initialization')
    endif
    close(input_file)
    ! debug
    ! write(*,*) trim(adjustl(mech_name))
    ! write(*,*) init_pressure
    ! write(*,*) init_temperature
    ! write(*,*) trim(adjustl(mass_or_mole))
    ! write(*,*) trim(adjustl(init_comp))

    mech_file = trim(mech_path)//'/'//trim(mech_name)//'.xml'
    ! load chemistry mechanism
    write(*,'(a)') '  read in the gas phase mechanism'
    write(*,'(a)') '  cantera file: '//trim(mech_name)//'.xml'
    gas = importphase(trim(mech_file), trim(mech_name)//'_mix')

    nspec = nspecies(gas)      ! number of species
    nreac = nreactions(gas)    ! number of reactions
    nelem = nelements(gas)     ! number of reactions
    ! ! debug
    ! write(*,*) 'number of elements is ', nelem
    ! write(*,*) 'number of species is ', nspec
    ! write(*,*) 'number of reactions is ', nreac

    call lisol_allocate()

    do i = 1, nelem
        call getElementName(gas, i, elem_name(i))
        ! debug
        ! write(*,*) 'Element ', i, ' is ', elem_name(i)
    enddo
    do i = 1, nspec
        call getSpeciesName(gas, i, spec_name(i))
        ! debug
        ! write(*,*) 'Species ', i, ' is ', spec_name(i)
    enddo
    do i = 1, nreac
        call getReactionString(gas, i, reac_name(i))
        ! reac_type(i) = reactionType(gas, i)
        itemp = isReversible(gas, i-1)  ! start from 0
        if(itemp .ne. 0)then
            isRevers(i) = .true.
        elseif(itemp .eq. 0)then
            isRevers(i) = .false.
        endif
        ! debug
        ! write(*,*) 'Reaction ', i, ' is ', reac_name(i)
    enddo
    ! ! debug
    ! write(*,*) 'Reaction type of reaction 7 is ', reac_type(7)
    do i = 1,nspec
        do j = 1,nreac
            reac_coeffs(i,j) = reactantStoichCoeff(gas, i, j)
            prod_coeffs(i,j) = productStoichCoeff(gas, i, j)
        enddo
    enddo
    ! ! debug
    ! write(*,*) 'Coefficient of species 1 in reactants among reaction 1 is ', reac_coeffs(1,1)
    ! write(*,*) 'Coefficient of species 1 in products among reaction 1 is ', prod_coeffs(1,1)
    ! write(*,*) 'Reactant coefficient in reaction 1 is', reac_coeffs(:,1)
    !*********************************************************************
    !> output these basic properties to files through Python interface 
    !> in Cantera, and then read them in formats here. We just need the
    !> data with the same order of precision in mechanism, right?

    write(*,'(a)') '  read in the properties of species'
    dir_comm = 'mkdir -p '//trim(inter_path)
	!--delete any existing directory
	dir_comm1 = 'rm -rf '//trim(inter_path)
	call system(trim(dir_comm1))
	!-and create a new one
	call system(trim(dir_comm))

    write(mech_file_length, '(i10)') len(trim(adjustl(mech_file)))
    open(interface_file,file=trim(fort2py),status='unknown',form='formatted', &
                        position='rewind',action='write')
    write(interface_file, '(a'//trim(adjustl(mech_file_length))//')') trim(adjustl(mech_file))
    close(interface_file)

    dir_comm = 'python '//trim(py_file)
    call execute_command_line(trim(dir_comm), wait = .true., exitstat = ierr)
    if(ierr == 0)then
        ! read the element atom number for each species
        do i = 1, nspec
            file_temp = trim(inter_path)//'/'//trim(spec_name(i))//spec_file_suffix
            open(interface_file,file=trim(file_temp),status='unknown',form='formatted', &
                                position='rewind',action='read')
            do j = 1, nelem
                read(interface_file, *) elem_num(i,j)
            enddo
            close(interface_file)
            ! ! debug
            ! if(spec_name(i) .eq. 'H2O')then
            !     do j = 1, nelem
            !         write(*, '(a2, i1)') elem_name(j), elem_num(i,j)
            !     enddo
            ! endif
        enddo
        ! read thermodynamic data for each species
        do i = 1, nspec
            file_temp = trim(inter_path)//'/'//trim(spec_name(i))//thermo_file_suffix
            open(interface_file,file=trim(file_temp),status='unknown',form='formatted', &
                                position='rewind',action='read')
            read(interface_file, *) spec_therm(i, 1)  ! low T
            read(interface_file, *) spec_therm(i, 2)  ! middle T
            read(interface_file, *) spec_therm(i, 3)  ! high T
            read(interface_file, *) spec_therm(i, 4:10)  ! coeffs in low T range
            read(interface_file, *) spec_therm(i, 11:spec_therm_num)  ! coeffs in high T range
            close(interface_file)
            ! ! debug
            ! if(spec_name(i) .eq. 'H2O')then
            !     write(*,*) 'Thermo of '//trim(spec_name(i))//':'
            !     write(*,*) spec_therm(i,1)
            !     write(*,*) spec_therm(i,2)
            !     write(*,*) spec_therm(i,3)
            !     write(*,'(7e14.8)') (spec_therm(i,j), j = 4, 10)
            !     write(*,'(7e14.8)') (spec_therm(i,j), j = 11, spec_therm_num)
            ! endif
        enddo
        ! read transport data for each species
        do i = 1, nspec
            file_temp = trim(inter_path)//'/'//trim(spec_name(i))//trans_file_suffix
            open(interface_file,file=trim(file_temp),status='unknown',form='formatted', &
                                position='rewind',action='read')
            do j = 1, spec_trans_num
                read(interface_file, *) spec_trans(i, j)
            enddo
            close(interface_file)
            ! ! debug
            ! if(spec_name(i) .eq. 'H2O')then
            !     write(*,*) 'Trans of '//trim(spec_name(i))//':'
            !     write(*,'(6f8.3)') (spec_trans(i,j), j = 1, spec_trans_num)
            ! endif
        enddo
        ! ISO units
        spec_trans_ISO(:,:) = spec_trans(:,:)
        do i = 1,nspec
            spec_trans_ISO(i,2) = spec_trans_ISO(i,2) * Boltzmann !< K --> J
            spec_trans_ISO(i,3) = spec_trans_ISO(i,3) / m2A  !< Angstroms --> m
            spec_trans_ISO(i,4) = Debye2CoulombM * spec_trans_ISO(i,4)  !< Debye --> Coulomb-m
            spec_trans_ISO(i,5) = spec_trans_ISO(i,5) / (m2A**3.d0)  !< A^3 --> m^3
        enddo
        !---- read in reaction types and corresponding parameters
        write(*,'(a)') '  read in the parameters of reactions'
        file_temp = trim(inter_path)//'/'//Total_reac_type
        open(interface_file,file=trim(file_temp),status='unknown',form='formatted', &
                            position='rewind',action='read')
        do i = 1,nreac
            read(interface_file,*) reac_type(i,1), reac_type(i,2), reac_type(i,3)
            ! write(*,*) reac_type(i,1), reac_type(i,2), reac_type(i,3)
        enddo
        ! get the max_Pnum from file reac_type.dat
        read(interface_file,*) max_Pnum
        if(max_Pnum .gt. 0) allocate(rate_coeff_PLOG(nreac, max_Pnum, 4))
        close(interface_file)
        !-- set the default values of 0
        rate_coeff(:,:) = 0.d0
        collision_efficiency(:,:) = 1.d0
        Pnum_PLOG(:) = 0
        if(max_Pnum .gt. 0) rate_coeff_PLOG(:,:,:) = 0.d0
        do i = 1,nreac
            write(ctemp, '(i10)') i
            ! write(*,*) 'reaction', i
            file_temp = trim(inter_path)//'/'//trim(adjustl(ctemp))//reac_file_suffix
            open(interface_file,file=trim(file_temp),status='unknown',form='formatted', &
                                position='rewind',action='read')
            ! 1. elementary reaction with Arrhenius rate
            if(reac_type(i,1) .eq. ELEMENTARY_RXN)then
                read(interface_file,*) rate_coeff(i,1), rate_coeff(i,2), rate_coeff(i,3) ! pre_exponential_factor, temperature_exponent, and activation_energy
            ! 2. three-body reaction with specific efficiency, else default 1.0
            elseif(reac_type(i,1) .eq. THREE_BODY_RXN)then
                read(interface_file,*) rate_coeff(i,1), rate_coeff(i,2), rate_coeff(i,3) ! pre_exponential_factor, temperature_exponent, and activation_energy
                if(reac_type(i,2) .ne. non)then
                    do j = 1,reac_type(i,2) ! number of collision efficiency parameters
                        read(interface_file,*) ctemp
                        itemp = speciesIndex(gas, trim(adjustl(ctemp)))
                        read(interface_file,*) collision_efficiency(i, itemp)
                    enddo
                endif
                ! ! debug
                ! if(i == 1)then
                !     write(*,*) collision_efficiency(i,:)
                ! endif
            ! 3. simple falloff reaction
            elseif((reac_type(i,1) .eq. FALLOFF_RXN) .and. (reac_type(i,2) .eq. SIMPLE_FALLOFF))then
                read(interface_file,*) rate_coeff(i,1), rate_coeff(i,2), rate_coeff(i,3) ! low pressure
                read(interface_file,*) rate_coeff(i,4), rate_coeff(i,5), rate_coeff(i,6) ! high pressure
                if(reac_type(i,3) .ne. non)then
                    do j = 1,reac_type(i,3) ! number of collision efficiency parameters
                        read(interface_file,*) ctemp
                        itemp = speciesIndex(gas, trim(adjustl(ctemp)))
                        read(interface_file,*) collision_efficiency(i, itemp)
                    enddo
                endif
            ! 4. TROE falloff reaction
            elseif((reac_type(i,1) .eq. FALLOFF_RXN) .and. (reac_type(i,2) .eq. TROE_FALLOFF))then
                read(interface_file,*) rate_coeff(i,1), rate_coeff(i,2), rate_coeff(i,3) ! low pressure
                read(interface_file,*) rate_coeff(i,4), rate_coeff(i,5), rate_coeff(i,6) ! high pressure
                read(interface_file,*) rate_coeff(i,7), rate_coeff(i,8), &
                                       rate_coeff(i,9), rate_coeff(i,10) ! TROE parameters
                if(reac_type(i,3) .ne. non)then
                    do j = 1,reac_type(i,3) ! number of collision efficiency parameters
                        read(interface_file,*) ctemp
                        itemp = speciesIndex(gas, trim(adjustl(ctemp)))
                        read(interface_file,*) collision_efficiency(i, itemp)
                    enddo
                endif
            ! 5. SRI-3 falloff reaction
            elseif((reac_type(i,1) .eq. FALLOFF_RXN) .and. (reac_type(i,2) .eq. SRI_FALLOFF_3))then
                read(interface_file,*) rate_coeff(i,1), rate_coeff(i,2), rate_coeff(i,3) ! low pressure
                read(interface_file,*) rate_coeff(i,4), rate_coeff(i,5), rate_coeff(i,6) ! high pressure
                read(interface_file,*) rate_coeff(i,7), rate_coeff(i,8), rate_coeff(i,9) ! SRI3 parameters
                rate_coeff(i,10) = 1.d0 ! d
                rate_coeff(i,11) = 0.d0 ! e
                if(reac_type(i,3) .ne. non)then
                    do j = 1,reac_type(i,3) ! number of collision efficiency parameters
                        read(interface_file,*) ctemp
                        itemp = speciesIndex(gas, trim(adjustl(ctemp)))
                        read(interface_file,*) collision_efficiency(i, itemp)
                    enddo
                endif
            ! 6. SRI-5 falloff reaction
            elseif((reac_type(i,1) .eq. FALLOFF_RXN) .and. (reac_type(i,2) .eq. SRI_FALLOFF_5))then
                read(interface_file,*) rate_coeff(i,1), rate_coeff(i,2), rate_coeff(i,3) ! low pressure
                read(interface_file,*) rate_coeff(i,4), rate_coeff(i,5), rate_coeff(i,6) ! high pressure
                read(interface_file,*) rate_coeff(i,7), rate_coeff(i,8), rate_coeff(i,9),&
                                       rate_coeff(i,10), rate_coeff(i,11) ! SRI5 parameters
                if(reac_type(i,3) .ne. non)then
                    do j = 1,reac_type(i,3) ! number of collision efficiency parameters
                        read(interface_file,*) ctemp
                        itemp = speciesIndex(gas, trim(adjustl(ctemp)))
                        read(interface_file,*) collision_efficiency(i, itemp)
                    enddo
                endif
            ! 7. PLOG reactions
            elseif(reac_type(i,1) .eq. PLOG_RXN)then
                do j = 1,reac_type(i,2) ! number of pressure points
                    Pnum_PLOG(i) = reac_type(i,2)
                    read(interface_file,*) rate_coeff_PLOG(i,j,1), rate_coeff_PLOG(i,j,2), &
                                           rate_coeff_PLOG(i,j,3), rate_coeff_PLOG(i,j,4) ! pressure, pre_exponential_factor, temperature_exponent, activation_energy
                enddo
            else
                call lisol_error(1, &
                'lisol: error reac_type, not supported yet! Only for elementary, Three-body, Falloff, and PLOG reactions now.')
            endif
            close(interface_file)
        enddo
        !----
    else
        call lisol_error(ierr, 'lisol: error during calling python !')
    endif

    ! finish the data exchange, and then delete the temporal path and file
    !--delete any existing directory
	dir_comm1 = 'rm -rf '//trim(inter_path)//' '//trim(fort2py)
    call system(trim(dir_comm1))
    !*********************************************************************

    !*********************************************************************
    !> set the solver variables, and test for subroutines
    ! ! call getAtomicWeights(gas, lisol_molmas)
    ! ! lisol_molmas(:) = lisol_molmas(:) * 1.d-3  !< convert kg/kmol to kg/mol
    ! ! write(*,*) 'Elements:'
    ! ! write(*,*) lisol_molmas
    call getMolecularWeights(gas, lisol_molmas) !< kg/kmol
    ! lisol_molmas(:) = lisol_molmas(:) * 1.d-3  !< convert kg/kmol to kg/mol
    ! ! debug
    ! ! write(*,*) 'Species:'
    ! ! write(*,*) lisol_molmas
    ! just test
    lisol_pressure(:,:,:) = init_pressure
    lisol_temperature(:,:,:) = init_temperature
    ! temp = (5000.d0 - 300.d0) / dble(dimx_total-1)
    ! do i = 0,dimx_total-1
    !     lisol_temperature(i+1,1,1) = 300.d0 + temp*i
    ! enddo

    init_frac(:) = 0.d0
    call split_string()
    ! normalized
    init_frac(:) = init_frac(:) / sum(init_frac(:))
    ! write(*,*) 'Init fraction is ', init_frac

    ! mass_or_mole = mole_description
    if(trim(mass_or_mole) .eq. mass_description)then
        lisol_mass_frac(:,:,:,:) = 0.d0
        do l = 1, nspec
            lisol_mass_frac(:,:,:,l) = init_frac(l)
        enddo
        call lisol_compute_meanmolmas('MassFrac')
        ! write(*,*) 'Mole fraction: ', lisol_mole_frac(1,1,1,:)
    elseif(trim(mass_or_mole) .eq. mole_description)then
        lisol_mole_frac(:,:,:,:) = 0.d0
        do l = 1, nspec
            lisol_mole_frac(:,:,:,l) = init_frac(l)
        enddo
        call lisol_compute_meanmolmas('MoleFrac')
        ! write(*,*) 'Mass fraction: ', lisol_mass_frac(1,1,1,:)
    else
        call lisol_error(1, 'lisol: error fraction type')
    endif
    
    ! debug
    ! write(*,*) 'Mean molecular weight is ', lisol_meanmolmas(1,1,1)
    ! write(*,*) 'Mean molecular weight is ', lisol_meanmolmas(2,1,1)
    call collision_integral_para_init()
    call setupCollisionIntegral()

    call lisol_compute_density()
    call lisol_compute_cp_species()
    call lisol_compute_cp()
    call lisol_compute_cv()
    call lisol_compute_enthalpy_species()
    call lisol_compute_enthalpy()
    call lisol_compute_gamma()
    call lisol_compute_lambda()
    call lisol_compute_viscosity()
    call lisol_compute_mixDiff()
    call lisol_compute_source_terms()
    !*********************************************************************

    call lisol_deallocate()
    ! stop
end program lisol
!=============================================
! This subroutine is used for allocate the
! array for the solver
!=============================================
subroutine lisol_allocate()
    use lisol_chemistry_mod
    implicit none
    integer(kind=4) :: i, s
    
    ! kinetic
    allocate(elem_name(nelem), spec_name(nspec))
    allocate(reac_name(nreac), reac_type(nreac, 3), reac_coeffs(nspec, nreac), prod_coeffs(nspec, nreac))
    allocate(isRevers(nreac), Pnum_PLOG(nreac))
    allocate(rate_coeff(nreac,11), collision_efficiency(nreac, nspec))
    allocate(lisol_heat_prod_rate(dimx_total, dimy_total, dimz_total))
    allocate(lisol_mass_prod_rate(dimx_total, dimy_total, dimz_total, nspec))
    allocate(lisol_heat_release(dimx_total, dimy_total, dimz_total))
    ! thermo
    allocate(elem_num(nspec, nelem), spec_therm(nspec, spec_therm_num), spec_trans(nspec, spec_trans_num))
    allocate(spec_trans_ISO(nspec, spec_trans_num))
    allocate(lisol_mass_frac(dimx_total, dimy_total, dimz_total, nspec))
    allocate(lisol_mole_frac(dimx_total, dimy_total, dimz_total, nspec))
    allocate(lisol_pressure(dimx_total, dimy_total, dimz_total))
    allocate(lisol_temperature(dimx_total, dimy_total, dimz_total))
    allocate(lisol_meanmolmas(dimx_total, dimy_total, dimz_total))
    allocate(lisol_density(dimx_total, dimy_total, dimz_total))
    allocate(lisol_cp_species(dimx_total, dimy_total, dimz_total, nspec))
    allocate(lisol_cp(dimx_total, dimy_total, dimz_total))
    allocate(lisol_cv(dimx_total, dimy_total, dimz_total))
    allocate(lisol_gamma(dimx_total, dimy_total, dimz_total))
    allocate(lisol_enthalpy(dimx_total, dimy_total, dimz_total))
    allocate(lisol_enthalpy_spec(dimx_total, dimy_total, dimz_total, nspec))
    allocate(lisol_molmas(nspec), init_frac(nspec))
    ! trans
    allocate(lisol_lambda(dimx_total, dimy_total, dimz_total))
    allocate(lisol_viscosity(dimx_total, dimy_total, dimz_total))
    allocate(m_epsilon(nspec, nspec), m_delta(nspec, nspec), m_reducedMass(nspec, nspec))
    m_epsilon(:,:) = 0.d0
    m_delta(:,:) = 0.d0
    m_reducedMass(:,:) = 0.d0
    allocate(m_diam(nspec, nspec), m_dipole(nspec, nspec))
    m_diam(:,:) = 0.d0
    allocate(m_crot(nspec), m_polar(nspec))
    allocate(m_visccoeffs(nspec, fitTrans_order+1), m_condcoeffs(nspec, fitTrans_order+1))
    s = 0
    do i = 1,nspec
        s = s + i
    enddo
    ! write(*,*) 's = ', s
    allocate(m_diffcoeffs(s, fitTrans_order+1))
    allocate(lisol_mixDiff(dimx_total, dimy_total, dimz_total, nspec))

    ! others
    allocate(m_bdiff(dimx_total, dimy_total, dimz_total, nspec, nspec))
    allocate(lisol_work(dimx_total, dimy_total, dimz_total), lisol_work2(dimx_total, dimy_total, dimz_total))
    allocate(lisol_work3(dimx_total, dimy_total, dimz_total), lisol_work4(dimx_total, dimy_total, dimz_total))
    allocate(lisol_work5(dimx_total, dimy_total, dimz_total))
    allocate(lisol_props_species(dimx_total, dimy_total, dimz_total, nspec))
    allocate(lisol_props_species2(dimx_total, dimy_total, dimz_total, nspec))
    allocate(flag(dimx_total, dimy_total, dimz_total))
    allocate(devia_coeffs(nspec))

    return
end subroutine lisol_allocate
!=============================================
! This subroutine is used for deallocate the
! array for the solver
!=============================================
subroutine lisol_deallocate()
    use lisol_chemistry_mod
    implicit none
    ! if(allocated()) deallocate()
    ! kinetic
    if(allocated(elem_name)) deallocate(elem_name)
    if(allocated(spec_name)) deallocate(spec_name)
    if(allocated(reac_name)) deallocate(reac_name)
    if(allocated(elem_num)) deallocate(elem_num)
    if(allocated(reac_type)) deallocate(reac_type)
    if(allocated(reac_coeffs)) deallocate(reac_coeffs)
    if(allocated(prod_coeffs)) deallocate(prod_coeffs)
    if(allocated(isRevers)) deallocate(isRevers)
    if(allocated(collision_efficiency)) deallocate(collision_efficiency)
    if(allocated(rate_coeff)) deallocate(rate_coeff)
    if(allocated(rate_coeff_PLOG)) deallocate(rate_coeff_PLOG)
    if(allocated(Pnum_PLOG)) deallocate(Pnum_PLOG)
    if(allocated(lisol_heat_prod_rate)) deallocate(lisol_heat_prod_rate)
    if(allocated(lisol_mass_prod_rate)) deallocate(lisol_mass_prod_rate)
    if(allocated(lisol_heat_release)) deallocate(lisol_heat_release)
    ! thermo
    if(allocated(spec_therm)) deallocate(spec_therm)
    if(allocated(lisol_mass_frac)) deallocate(lisol_mass_frac)
    if(allocated(lisol_mole_frac)) deallocate(lisol_mole_frac)
    if(allocated(lisol_pressure)) deallocate(lisol_pressure)
    if(allocated(lisol_temperature)) deallocate(lisol_temperature)
    if(allocated(lisol_meanmolmas)) deallocate(lisol_meanmolmas)
    if(allocated(lisol_density)) deallocate(lisol_density)
    if(allocated(lisol_cp_species)) deallocate(lisol_cp_species)
    if(allocated(lisol_cp)) deallocate(lisol_cp)
    if(allocated(lisol_cv)) deallocate(lisol_cv)
    if(allocated(lisol_gamma)) deallocate(lisol_gamma)
    if(allocated(lisol_enthalpy)) deallocate(lisol_enthalpy)
    if(allocated(lisol_enthalpy_spec)) deallocate(lisol_enthalpy_spec)
    if(allocated(lisol_molmas)) deallocate(lisol_molmas)
    if(allocated(init_frac)) deallocate(init_frac)
    ! trans
    if(allocated(spec_trans)) deallocate(spec_trans)
    if(allocated(spec_trans_ISO)) deallocate(spec_trans_ISO)
    if(allocated(lisol_lambda)) deallocate(lisol_lambda)
    if(allocated(lisol_viscosity)) deallocate(lisol_viscosity)
    if(allocated(m_epsilon)) deallocate(m_epsilon)
    if(allocated(m_delta)) deallocate(m_delta)
    if(allocated(m_reducedMass)) deallocate(m_reducedMass)
    if(allocated(m_diam)) deallocate(m_diam)
    if(allocated(m_crot)) deallocate(m_crot)
    if(allocated(m_dipole)) deallocate(m_dipole)
    if(allocated(m_polar)) deallocate(m_polar)
    if(allocated(m_visccoeffs)) deallocate(m_visccoeffs)
    if(allocated(m_condcoeffs)) deallocate(m_condcoeffs)
    if(allocated(m_diffcoeffs)) deallocate(m_diffcoeffs)
    if(allocated(lisol_mixDiff)) deallocate(lisol_mixDiff)
    ! others
    if(allocated(flag)) deallocate(flag)
    if(allocated(lisol_work)) deallocate(lisol_work)
    if(allocated(lisol_work2)) deallocate(lisol_work2)
    if(allocated(lisol_work3)) deallocate(lisol_work3)
    if(allocated(lisol_work4)) deallocate(lisol_work4)
    if(allocated(lisol_work5)) deallocate(lisol_work5)
    if(allocated(lisol_props_species)) deallocate(lisol_props_species)
    if(allocated(lisol_props_species2)) deallocate(lisol_props_species2)
    if(allocated(m_bdiff)) deallocate(m_bdiff)
    if(allocated(devia_coeffs)) deallocate(devia_coeffs)

    return
end subroutine lisol_deallocate
!=============================================
! This subroutine is used for split the string
! array of species
!=============================================
subroutine split_string()
    use lisol_chemistry_mod
    implicit none
    character(len=500) :: str
    character(len=50) :: delimiter
    character(len=200), allocatable, dimension(:) :: result
    integer(kind=4) :: count, i, j
    logical :: flag_spec
    
    integer(kind=4) :: pos, start, len_delim
    character(len=200) :: temp
    write(*,'(a)') '  Read in the initial fraction'
    allocate(result(nspec))

    ! split by ','
    str = trim(init_comp)
    delimiter = ','
    start = 1
    len_delim = len(trim(delimiter))
    count = 0
    do
        pos = index(str(start:), trim(delimiter))
        ! write(*,*) 'Pos: ', pos
        if (pos .eq. 0) exit

        temp = str(start:start+pos-2)
        count = count + 1
        result(count) = trim(temp)

        start = start + pos + len_delim - 1
    enddo
    ! write(*,*) 'start: ', start

    if(start .le. len(str)) then
        count = count + 1
        result(count) = trim(str(start:))
    endif

    ! write(*,*) 'count: ', count
    ! split by ':'
    delimiter = ':'
    do i = 1,count
        ! write(*,*) 'i = ', i, ', result = ', trim(result(i))
        temp = trim(result(i))
        pos = index(temp, trim(delimiter))
        flag_spec = .false.
        do j = 1,nspec
            if(trim(adjustl(temp(1:pos-1))) .eq. spec_name(j))then
                read(temp(pos+1:),*) init_frac(j)
                flag_spec = .true.
                ! write(*,*) 'String of init_frac is ', temp(pos+1:)
            endif
        enddo
        if(.not. flag_spec)then
            call lisol_error(1, 'split_string: error species name, not found in the mechanism')
        endif
    enddo

    deallocate(result)
    return
end subroutine split_string
!=============================================
! This subroutine is used to calculate the
! mean molarmass, and update mole fraction
! i.e., lisol_meanmolmas and lisol_mole_frac
!=============================================
subroutine lisol_compute_meanmolmas(description)
    use lisol_chemistry_mod
    implicit none
    character(len=*), intent(in) :: description
    integer(kind=4) :: i
    write(*,'(a)') '  Compute the mean molecular mass'
    if(trim(description) .eq. 'MassFrac')then
        do i = 1,nspec
            lisol_props_species(:,:,:,i) = lisol_mass_frac(:,:,:,i) / lisol_molmas(i)
        enddo
        lisol_meanmolmas(:,:,:) = sum(lisol_props_species(:,:,:,:), dim=4)
        lisol_meanmolmas(:,:,:) = 1.d0 / lisol_meanmolmas(:,:,:)
        do i = 1,nspec
            lisol_mole_frac(:,:,:,i) = lisol_props_species(:,:,:,i) * lisol_meanmolmas(:,:,:)
        enddo
    elseif(trim(description) .eq. 'MoleFrac')then
        do i = 1,nspec
            lisol_props_species(:,:,:,i) = lisol_mole_frac(:,:,:,i) * lisol_molmas(i)
        enddo
        lisol_meanmolmas(:,:,:) = sum(lisol_props_species(:,:,:,:), dim=4)
        do i = 1,nspec
            lisol_mass_frac(:,:,:,i) = lisol_props_species(:,:,:,i) / lisol_meanmolmas(:,:,:)
        enddo
    endif
    ! ! debug
    ! write(*,*) 'Mean molecular mass: ', lisol_meanmolmas(1,1,1)
    return
end subroutine lisol_compute_meanmolmas
!=============================================
! This subroutine is used for calculate the
! density
! i.e., lisol_density
!=============================================
subroutine lisol_compute_density()
    use lisol_chemistry_mod
    implicit none
    integer(kind=4) :: fid, i
    if(ideal_gas)then  ! ideal gas assumption
        write(*,'(a)') '  Compute the density with ideal-gas assumption'
        lisol_density(:,:,:) = lisol_pressure(:,:,:) * lisol_meanmolmas(:,:,:) &
                               / GasConstant / lisol_temperature(:,:,:)
    else
        call lisol_error(1, 'lisol_compute_density: not coded yet')
    endif

    ! open(fid,file=trim('temper_rho.dat'),status='unknown',form='formatted', &
    !                     position='rewind',action='write')
    ! write(fid, '(a)') trim('temper[K]          rho[kg/m3]')
    ! do i = 1,dimx_total
    !     write(fid, '(2e24.16)') lisol_temperature(i,1,1), lisol_density(i,1,1)
    ! enddo
    ! close(fid)
    ! debug
    ! write(*,*) 'Density: ', lisol_density(1,1,1)
    ! write(*,*) 'Density: ', lisol_density(dimx_total,dimy_total,dimz_total)
    return
end subroutine lisol_compute_density
!=============================================
! This subroutine is used for calculate the
! specific heat capacity at constant pressure
! for each species
! i.e., lisol_cp_species
!=============================================
subroutine lisol_compute_cp_species()
    use lisol_chemistry_mod
    implicit none
    integer(kind=4) :: l

    if(ideal_gas)then  ! ideal gas assumption
        write(*,'(a)') '  Compute the species specific heat capacity at constant pressure with ideal-gas assumption'
        lisol_work(:,:,:)  = (lisol_temperature(:,:,:)**2.d0)
        lisol_work2(:,:,:) = (lisol_temperature(:,:,:)**3.d0)
        lisol_work3(:,:,:) = (lisol_temperature(:,:,:)**4.d0)
        do l = 1,nspec
            flag = (lisol_temperature(:,:,:) <= spec_therm(l,2))  ! middle temperature
            lisol_cp_species(:,:,:,l) = GasConstant * &
            merge((spec_therm(l, 4) + spec_therm(l, 5)*lisol_temperature(:,:,:) &
                                    + spec_therm(l, 6)*lisol_work(:,:,:) &
                                    + spec_therm(l, 7)*lisol_work2(:,:,:) &
                                    + spec_therm(l, 8)*lisol_work3(:,:,:)), &
                  (spec_therm(l,11) + spec_therm(l,12)*lisol_temperature(:,:,:) &
                                    + spec_therm(l,13)*lisol_work(:,:,:) &
                                    + spec_therm(l,14)*lisol_work2(:,:,:) &
                                    + spec_therm(l,15)*lisol_work3(:,:,:)), flag)
            ! J/kmol/K / (kg/kmol) --> J/kg/K
            lisol_cp_species(:,:,:,l) = lisol_cp_species(:,:,:,l) / lisol_molmas(l)
            ! debug
            if(spec_name(l) == 'H2')then
                write(*,*) 'Cp_species of '//trim(spec_name(l))//' is ', lisol_cp_species(1,1,1,l)
                ! write(*,*) 'Cp_species of '//trim(spec_name(l))//' is ', lisol_cp_species(dimx_total, dimy_total, dimz_total, l)
            endif
        enddo
    else
        call lisol_error(1, 'lisol_compute_cp_species: not coded yet')
    endif
    ! write(*,*) 'Cp_species: ', lisol_cp_species(1,1,1,1)
    return
end subroutine lisol_compute_cp_species
!=============================================
! This subroutine is used for calculate the
! specific heat capacity at constant pressure
! i.e., lisol_cp, obtained from mass-weighted
!=============================================
subroutine lisol_compute_cp()
    use lisol_chemistry_mod
    implicit none
    integer(kind=4) :: l, i, fid

    if(ideal_gas)then  ! ideal gas assumption
        write(*,'(a)') '  Compute the specific heat capacity at constant pressure with ideal-gas assumption'
        lisol_cp(:,:,:) = sum(lisol_mass_frac(:,:,:,:) * lisol_cp_species(:,:,:,:), dim = 4)
        ! do l = 1,nspec
        !     lisol_cp(:,:,:) = lisol_cp(:,:,:) + lisol_mass_frac(:,:,:,l) * lisol_cp_species(:,:,:,l)
        ! enddo
        ! debug
        write(*,*) 'Cp: ', lisol_cp(1,1,1)
        ! write(*,*) 'Cp: ', lisol_cp(dimx_total, dimy_total, dimz_total)
    else
        call lisol_error(1, 'lisol_compute_cp: not coded yet')
    endif

    ! open(fid,file=trim('temper_cp.dat'),status='unknown',form='formatted', &
    !                     position='rewind',action='write')
    ! write(fid, '(a)') trim('temper[K]          cp[J/kg/K]')
    ! do i = 1,dimx_total
    !     write(fid, '(2e24.16)') lisol_temperature(i,1,1), lisol_cp(i,1,1)
    ! enddo
    ! close(fid)
    ! write(*,*) 'Cp: ', lisol_cp(1,1,1)
    return
end subroutine lisol_compute_cp
!=============================================
! This subroutine is used for calculate the
! specific heat capacity at constant volume
! i.e., lisol_cv
!=============================================
subroutine lisol_compute_cv()
    use lisol_chemistry_mod
    implicit none
    if(ideal_gas)then  ! ideal gas assumption
        write(*,'(a)') '  Compute the specific heat capacity at constant volume with ideal-gas assumption'
        lisol_cv(:,:,:) = lisol_cp(:,:,:) - GasConstant / lisol_meanmolmas(:,:,:)
        ! debug
        write(*,*) 'Cv: ', lisol_cv(1,1,1)
        ! write(*,*) 'Cv: ', lisol_cv(dimx_total, dimy_total, dimz_total)
    else
        call lisol_error(1, 'lisol_compute_cv: not coded yet')
    endif
    ! write(*,*) 'Cv: ', lisol_cv(1,1,1)
    return
end subroutine lisol_compute_cv
!=============================================
! This subroutine is used for calculate the
! enthalpy for each species
! i.e., lisol_enthalpy_spec
!=============================================
subroutine lisol_compute_enthalpy_species()
    use lisol_chemistry_mod
    implicit none
    integer(kind=4) :: l

    if(ideal_gas)then  ! ideal gas assumption
        write(*,'(a)') '  Compute the species enthalpy with ideal-gas assumption'
        lisol_work(:,:,:)  = lisol_temperature(:,:,:)/2.d0
        lisol_work2(:,:,:) = (lisol_temperature(:,:,:)**2.d0)/3.d0
        lisol_work3(:,:,:) = (lisol_temperature(:,:,:)**3.d0)/4.d0
        lisol_work4(:,:,:) = (lisol_temperature(:,:,:)**4.d0)/5.d0
        do l = 1,nspec
            flag = (lisol_temperature(:,:,:) <= spec_therm(l,2))  ! middle temperature
            lisol_enthalpy_spec(:,:,:,l) = GasConstant * lisol_temperature(:,:,:) *&
            merge((spec_therm(l, 4) + spec_therm(l, 5)*lisol_work(:,:,:) &
                                    + spec_therm(l, 6)*lisol_work2(:,:,:) &
                                    + spec_therm(l, 7)*lisol_work3(:,:,:) &
                                    + spec_therm(l, 8)*lisol_work4(:,:,:) &
                                    + spec_therm(l, 9)/lisol_temperature(:,:,:)), &
                  (spec_therm(l,11) + spec_therm(l,12)*lisol_work(:,:,:) &
                                    + spec_therm(l,13)*lisol_work2(:,:,:) &
                                    + spec_therm(l,14)*lisol_work3(:,:,:) &
                                    + spec_therm(l,15)*lisol_work4(:,:,:) &
                                    + spec_therm(l,16)/lisol_temperature(:,:,:)), flag)
            ! J/kmol / (kg/kmol) --> J/kg
            lisol_enthalpy_spec(:,:,:,l) = lisol_enthalpy_spec(:,:,:,l) / lisol_molmas(l)
            ! debug
            if(spec_name(l) == 'H2')then
                write(*,*) 'Enthalpy_species of '//trim(spec_name(l))//' is ', lisol_enthalpy_spec(1,1,1,l)
                ! write(*,*) 'Enthalpy_speciesof '//trim(spec_name(l))//' is ', lisol_enthalpy_spec(dimx_total, dimy_total, dimz_total, l)
            endif
        enddo
    else
        call lisol_error(1, 'lisol_compute_enthalpy_species: not coded yet')
    endif
    ! write(*,*) 'Enthalpy_species: ', lisol_enthalpy_spec(1,1,1,1)

    return
end subroutine lisol_compute_enthalpy_species
!=============================================
! This subroutine is used for calculate the
! enthalpy
! i.e., lisol_enthalpy
!=============================================
subroutine lisol_compute_enthalpy()
    use lisol_chemistry_mod
    implicit none
    integer(kind=4) :: l, i, fid

    if(ideal_gas)then  ! ideal gas assumption
        write(*,'(a)') '  Compute the enthalpy with ideal-gas assumption'
        lisol_enthalpy(:,:,:) = sum(lisol_mass_frac(:,:,:,:) * lisol_enthalpy_spec(:,:,:,:), dim = 4)
        ! do l = 1,nspec
        !     lisol_enthalpy(:,:,:) = lisol_enthalpy(:,:,:) + lisol_mass_frac(:,:,:,l) * lisol_enthalpy_spec(:,:,:,l)
        ! enddo
        ! debug
        write(*,*) 'Enthalpy: ', lisol_enthalpy(1,1,1)
        ! write(*,*) 'Enthalpy: ', lisol_enthalpy(dimx_total, dimy_total, dimz_total)
    else
        call lisol_error(1, 'lisol_compute_enthalpy: not coded yet')
    endif

    ! open(fid,file=trim('temper_enthalpy.dat'),status='unknown',form='formatted', &
    !                     position='rewind',action='write')
    ! write(fid, '(a)') trim('temper[K]          h[J/kg]')
    ! do i = 1,dimx_total
    !     write(fid, '(2e24.16)') lisol_temperature(i,1,1), lisol_enthalpy(i,1,1)
    ! enddo
    ! close(fid)

    return
end subroutine lisol_compute_enthalpy
!=============================================
! This subroutine is used for calculate the
! specific heat
! i.e., lisol_gamma
!=============================================
subroutine lisol_compute_gamma()
    use lisol_chemistry_mod
    implicit none
    write(*,'(a)') '  Compute the specific heat with ideal-gas assumption'
    lisol_gamma(:,:,:) = lisol_cp(:,:,:) / lisol_cv(:,:,:)
    ! ! debug
    ! write(*,*) 'Gamma: ', lisol_gamma(1,1,1)
    ! write(*,*) 'Gamma: ', lisol_gamma(dimx_total, dimy_total, dimz_total)
    return
end subroutine lisol_compute_gamma
!=============================================
! This subroutine is used for calculate the
! thermal conductivity
! i.e., lisol_lambda
!=============================================
subroutine lisol_compute_lambda()
    use lisol_chemistry_mod
    implicit none
    integer(kind=4) :: l, i, fid
    write(*,'(a)') '  Compute the thermal conductivity with ideal-gas assumption'
    ! same as Cantera, using the fitted 4th order formula, for species k, there is
    ! lambda(k) = T^{1/2} * \sum_{n=0}^{4}{b_n(k) (\ln T)^n}, where b_n(k) are the coeffs for species k
    lisol_work(:,:,:)  = dsqrt(lisol_temperature(:,:,:))
    lisol_work2(:,:,:) = dlog(lisol_temperature(:,:,:))
    lisol_work3(:,:,:) = (dlog(lisol_temperature(:,:,:)))**(2.d0)
    lisol_work4(:,:,:) = (dlog(lisol_temperature(:,:,:)))**(3.d0)
    lisol_work5(:,:,:) = (dlog(lisol_temperature(:,:,:)))**(4.d0)
    do l = 1,nspec
        lisol_props_species(:,:,:,l) = lisol_work(:,:,:) * ( m_condcoeffs(l,1) &
                                                           + m_condcoeffs(l,2)*lisol_work2(:,:,:) &
                                                           + m_condcoeffs(l,3)*lisol_work3(:,:,:) &
                                                           + m_condcoeffs(l,4)*lisol_work4(:,:,:) &
                                                           + m_condcoeffs(l,5)*lisol_work5(:,:,:) )
        ! lisol_work(:,:,:) = lisol_temperature(:,:,:)**(0.5d0) * ( m_condcoeffs(l,1) &
        !                     + m_condcoeffs(l,2)*(dlog(lisol_temperature(:,:,:))) &
        !                     + m_condcoeffs(l,3)*((dlog(lisol_temperature(:,:,:)))**(2.d0)) &
        !                     + m_condcoeffs(l,4)*((dlog(lisol_temperature(:,:,:)))**(3.d0)) &
        !                     + m_condcoeffs(l,5)*((dlog(lisol_temperature(:,:,:)))**(4.d0)) )
        ! lisol_props_species(:,:,:,l) = lisol_work(:,:,:)
        ! debug
        if(spec_name(l) == 'H2')then
            write(*,*) 'Thermal conductivity of '//trim(spec_name(l))//' is ', lisol_props_species(1,1,1,l)
            ! write(*,*) 'Thermal conductivity of '//trim(spec_name(l))//' is ', lisol_props_species(dimx_total, dimy_total, dimz_total, l)
        endif
    enddo

    ! now calculate the mixture-averaged thermal conductivity
    !!!!! NOTE, the mole fraction needs to be updated
    lisol_work(:,:,:)  = sum(lisol_mole_frac(:,:,:,:) * lisol_props_species(:,:,:,:), dim = 4)
    lisol_work2(:,:,:) = sum(lisol_mole_frac(:,:,:,:) / lisol_props_species(:,:,:,:), dim = 4)
    ! do l = 1,nspec
    !     lisol_work(:,:,:) = lisol_work(:,:,:) + lisol_mole_frac(:,:,:,l) * lisol_props_species(:,:,:,l)
    !     lisol_work2(:,:,:) = lisol_work2(:,:,:) + lisol_mole_frac(:,:,:,l) / lisol_props_species(:,:,:,l)
    ! enddo
    lisol_lambda(:,:,:) = 0.5d0 * (lisol_work(:,:,:) + 1.d0 / lisol_work2(:,:,:))
    ! debug
    write(*,*) 'Thermal conductivity is ', lisol_lambda(1,1,1)
    ! write(*,*) 'Thermal conductivity is ', lisol_lambda(dimx_total, dimy_total, dimz_total)

    ! open(fid,file=trim('temper_lambda.dat'),status='unknown',form='formatted', &
    !                     position='rewind',action='write')
    ! write(fid, '(a)') trim('temper[K]          lambda[W/m/K]')
    ! do i = 1,dimx_total
    !     write(fid, '(2e24.16)') lisol_temperature(i,1,1), lisol_lambda(i,1,1)
    ! enddo
    ! close(fid)

    return
end subroutine lisol_compute_lambda
!=============================================
! This subroutine is used for calculate the
! viscosity
! i.e., lisol_viscosity
!=============================================
subroutine lisol_compute_viscosity()
    use lisol_chemistry_mod
    implicit none
    integer(kind=4) :: l, i, k, j, fid
    write(*,'(a)') '  Compute the viscosity with ideal-gas assumption'
    ! same as Cantera, using the fitted 4th order formula, for species k, there is
    ! (eta(k))^{1/2} = T^{1/4} * \sum_{n=0}^{4}{a_n(k) (\ln T)^n}, where a_n(k) are the coeffs for species k
    lisol_work(:,:,:)  = lisol_temperature(:,:,:)**(0.25d0)
    lisol_work2(:,:,:) = dlog(lisol_temperature(:,:,:))
    lisol_work3(:,:,:) = (dlog(lisol_temperature(:,:,:)))**(2.d0)
    lisol_work4(:,:,:) = (dlog(lisol_temperature(:,:,:)))**(3.d0)
    lisol_work5(:,:,:) = (dlog(lisol_temperature(:,:,:)))**(4.d0)
    do l = 1,nspec
        lisol_props_species(:,:,:,l) = (lisol_work(:,:,:) * ( m_visccoeffs(l,1) &
                                                            + m_visccoeffs(l,2)*lisol_work2(:,:,:) &
                                                            + m_visccoeffs(l,3)*lisol_work3(:,:,:) &
                                                            + m_visccoeffs(l,4)*lisol_work4(:,:,:) &
                                                            + m_visccoeffs(l,5)*lisol_work5(:,:,:) ))**(2.d0)
        ! lisol_work(:,:,:) = lisol_temperature(:,:,:)**(0.25d0) * ( m_visccoeffs(l,1) &
        !                     + m_visccoeffs(l,2)*(dlog(lisol_temperature(:,:,:))) &
        !                     + m_visccoeffs(l,3)*((dlog(lisol_temperature(:,:,:)))**(2.d0)) &
        !                     + m_visccoeffs(l,4)*((dlog(lisol_temperature(:,:,:)))**(3.d0)) &
        !                     + m_visccoeffs(l,5)*((dlog(lisol_temperature(:,:,:)))**(4.d0)) )
        ! lisol_props_species(:,:,:,l) = (lisol_work(:,:,:))**(2.d0)
        ! debug
        if(trim(spec_name(l)) == 'H2')then
            write(*,*) 'Viscosity of '//trim(spec_name(l))//' is ', lisol_props_species(1,1,1,l)
            ! write(*,*) 'Viscosity of '//trim(spec_name(l))//' is ', lisol_props_species(dimx_total, dimy_total, dimz_total, l)
        endif
    enddo
    ! now calculate the mixture-averaged viscosity
    !!!!! NOTE, the mole fraction needs to be updated
    lisol_viscosity(:,:,:) = 0.d0
    do k = 1,nspec
        lisol_work(:,:,:) = 0.d0
        do j = 1,nspec
            lisol_work2(:,:,:) = 1.d0/dsqrt(8.d0) * (1.d0 + lisol_molmas(k)/lisol_molmas(j))**(-0.5d0) &
                                 * (1.d0 + dsqrt(lisol_props_species(:,:,:,k)/lisol_props_species(:,:,:,j)) &
                                           * (lisol_molmas(j)/lisol_molmas(k))**(1.d0/4.d0) )**(2.d0) !< PHI_kj
            lisol_work(:,:,:) = lisol_work(:,:,:) + lisol_mole_frac(:,:,:,j) * lisol_work2(:,:,:)
        enddo
        lisol_viscosity(:,:,:) = lisol_viscosity(:,:,:) + lisol_mole_frac(:,:,:,k) &
                                                        * lisol_props_species(:,:,:,k) / lisol_work(:,:,:)
    enddo
    ! debug
    write(*,*) 'Viscosity is ', lisol_viscosity(1,1,1)
    ! write(*,*) 'Viscosity is ', lisol_viscosity(dimx_total, dimy_total, dimz_total)

    ! open(fid,file=trim('temper_viscosity.dat'),status='unknown',form='formatted', &
    !                     position='rewind',action='write')
    ! write(fid, '(a)') trim('temper[K]          mu[Pa*s]')
    ! do i = 1,dimx_total
    !     write(fid, '(2e24.16)') lisol_temperature(i,1,1), lisol_viscosity(i,1,1)
    ! enddo
    ! close(fid)

    return
end subroutine lisol_compute_viscosity
!=============================================
! This subroutine is used for calculate the
! mixture-averaged diffusion coefficient
! i.e., lisol_mixDiff
! NOTE!!!!! correction velocity should be considered later
!=============================================
subroutine lisol_compute_mixDiff()
    use lisol_chemistry_mod
    implicit none
    integer(kind=4) :: icount, i, j, k, fid
    character(len=10) :: temp
    write(*,'(a)') '  Compute the mixture-averaged diffusion coefficients'
    ! 1 --> calculate the binary diffusion coefficients at unit pressure
    lisol_work(:,:,:)  = dsqrt(lisol_temperature(:,:,:))
    lisol_work2(:,:,:) = dlog(lisol_temperature(:,:,:))
    lisol_work3(:,:,:) = (dlog(lisol_temperature(:,:,:)))**(2.d0)
    lisol_work4(:,:,:) = (dlog(lisol_temperature(:,:,:)))**(3.d0)
    lisol_work5(:,:,:) = (dlog(lisol_temperature(:,:,:)))**(4.d0)
    icount = 1
    do i = 1,nspec
        do j = i,nspec
            m_bdiff(:,:,:,i,j) = lisol_temperature(:,:,:) * lisol_work(:,:,:) &
                                  * ( m_diffcoeffs(icount,1) &
                                    + m_diffcoeffs(icount,2)*lisol_work2(:,:,:) &
                                    + m_diffcoeffs(icount,3)*lisol_work3(:,:,:) &
                                    + m_diffcoeffs(icount,4)*lisol_work4(:,:,:) &
                                    + m_diffcoeffs(icount,5)*lisol_work5(:,:,:) )
            ! m_bdiff(:,:,:,i,j) = lisol_temperature(:,:,:) * dsqrt(lisol_temperature(:,:,:)) &
            !                       * ( m_diffcoeffs(icount,1) &
            !                         + m_diffcoeffs(icount,2)*dlog(lisol_temperature(:,:,:)) &
            !                         + m_diffcoeffs(icount,3)*(dlog(lisol_temperature(:,:,:)))**(2.d0) &
            !                         + m_diffcoeffs(icount,4)*(dlog(lisol_temperature(:,:,:)))**(3.d0) &
            !                         + m_diffcoeffs(icount,5)*(dlog(lisol_temperature(:,:,:)))**(4.d0) )
            m_bdiff(:,:,:,j,i) = m_bdiff(:,:,:,i,j)
            icount = icount + 1
        enddo
    enddo
    ! 2 --> get the mixture-averaged diffusion coefficients
    if(nspec .eq. 1)then
        lisol_mixDiff(:,:,:,1) = m_bdiff(:,:,:,1,1) / lisol_pressure(:,:,:)
    else
        do k = 1,nspec
            ! directly as the formula
            lisol_work(:,:,:) = 0.d0
            do j = 1,nspec
                if(j .ne. k)then
                    lisol_work(:,:,:) = lisol_work(:,:,:) + lisol_mole_frac(:,:,:,j)/m_bdiff(:,:,:,j,k)
                endif
            enddo
            ! where(lisol_work <= 0.d0)
            !     lisol_mixDiff(:,:,:,k) = m_bdiff(:,:,:,k,k) / lisol_pressure(:,:,:)
            ! elsewhere
            !     lisol_mixDiff(:,:,:,k) = (lisol_meanmolmas(:,:,:) - lisol_mole_frac(:,:,:,k)*lisol_molmas(k)) &
            !                              / (lisol_pressure(:,:,:) * lisol_meanmolmas(:,:,:) * lisol_work(:,:,:))
            ! end where
            lisol_mixDiff(:,:,:,k) = (1.d0 - lisol_mass_frac(:,:,:,k)) / lisol_work(:,:,:) / lisol_pressure(:,:,:)
            ! debug
            if(trim(spec_name(k)) == 'H2')then
                write(*,*) 'Mixture-averaged diffusion coeffs for '//trim(spec_name(k))//' is ', lisol_mixDiff(1,1,1,k)
                ! write(*,*) 'Mixture-averaged diffusion coeffs for '//trim(spec_name(k))//' is ', &
                !                                                     lisol_mixDiff(dimx_total, dimy_total, dimz_total,k)
            endif
        enddo
    endif

    ! write(temp, '(i5)') nspec+1
    ! open(fid,file=trim('temper_mixDiff.dat'),status='unknown',form='formatted', &
    !                     position='rewind',action='write')
    ! write(fid, '(a)') trim('temper[K]          Dk[m2/s] in the order of species index')
    ! do i = 1,dimx_total
    !     write(fid, '('//trim(adjustl(temp))//'e24.16)') lisol_temperature(i,1,1), (lisol_mixDiff(i,1,1,j), j = 1, nspec)
    ! enddo
    ! close(fid)

    return
end subroutine lisol_compute_mixDiff
!=============================================
! This subroutine is used for calculating the
! progress rates of each reaction firstly, then 
! the mass source terms for each species
!=============================================
subroutine lisol_compute_source_terms()
    use lisol_chemistry_mod
    implicit none
    integer(kind=4) :: i, j, itemp
    real(kind=8) :: rtemp

    write(*,'(a)') '  Compute the chemical potential and concentration for each species'
    lisol_work2(:,:,:) = dlog(lisol_temperature(:,:,:))
    lisol_work3(:,:,:) = (lisol_temperature(:,:,:)**2.d0)/2.d0
    lisol_work4(:,:,:) = (lisol_temperature(:,:,:)**3.d0)/3.d0
    lisol_work5(:,:,:) = (lisol_temperature(:,:,:)**4.d0)/4.d0
    do i = 1,nspec
        !---- get the species chemical potential
        ! firstly the entropy, J/kmol/K
        flag = (lisol_temperature(:,:,:) <= spec_therm(i,2))  ! middle temperature
        lisol_work(:,:,:) = GasConstant *&
          merge((spec_therm(i, 4)*lisol_work2(:,:,:) + spec_therm(i, 5)*lisol_temperature(:,:,:) &
                                                     + spec_therm(i, 6)*lisol_work3(:,:,:) &
                                                     + spec_therm(i, 7)*lisol_work4(:,:,:) &
                                                     + spec_therm(i, 8)*lisol_work5(:,:,:) &
                                                     + spec_therm(i,10)), &
                (spec_therm(i,11)*lisol_work2(:,:,:) + spec_therm(i,12)*lisol_temperature(:,:,:) &
                                                     + spec_therm(i,13)*lisol_work3(:,:,:) &
                                                     + spec_therm(i,14)*lisol_work4(:,:,:) &
                                                     + spec_therm(i,15)*lisol_work5(:,:,:) &
                                                     + spec_therm(i,17)), flag)
        lisol_props_species(:,:,:,i) = lisol_enthalpy_spec(:,:,:,i) * lisol_molmas(i) ! enthalpy, J/kmol
        lisol_props_species(:,:,:,i) = lisol_props_species(:,:,:,i) - lisol_work(:,:,:) * lisol_temperature(:,:,:) ! mu0 = G0 = H-TS
        !---- get the concentration, Ci = Xi*P/R/T, kmol/m3
        lisol_props_species2(:,:,:,i) = lisol_mole_frac(:,:,:,i) * lisol_pressure(:,:,:) &
                                        / GasConstant / lisol_temperature(:,:,:)
        ! ! debug
        ! if(spec_name(i) == 'CO')then
        !     write(*,*) 'Entropy [J/kmol/K] of ', trim(adjustl(spec_name(i))), ' is', lisol_work(1,1,1)
        !     write(*,*) 'Chemical potential [J/kmol] of ', trim(adjustl(spec_name(i))), ' is', lisol_props_species(1,1,1,i)
        !     write(*,*) 'Concentration [kmol/m3] of ', trim(adjustl(spec_name(i))), ' is', lisol_props_species2(1,1,1,i)
        ! endif
    enddo
    lisol_mass_prod_rate(:,:,:,:) = 0.d0
    lisol_heat_release(:,:,:) = dlog(OneAtm / GasConstant / lisol_temperature(:,:,:)) ! used as temp array here
    !---- calculate the rates of progress
    write(*,'(a)') '  Compute the mass production rates'
    do i = 1,nreac
        ! write(*,*) 'reaction', i
        devia_coeffs(:) = prod_coeffs(:,i) - reac_coeffs(:,i) ! v'' - v'
        !---- forward (or back) rate of progress
        ! 1. elementary reaction with Arrhenius rate
        if(reac_type(i,1) .eq. ELEMENTARY_RXN)then
            ! forward rate constant
            lisol_work(:,:,:) = rate_coeff(i,1) * lisol_temperature(:,:,:)**(rate_coeff(i,2)) &
                                * dexp(-rate_coeff(i,3) / GasConstant / lisol_temperature(:,:,:))
            lisol_work3(:,:,:) = lisol_work(:,:,:)
            ! debug
            if(if_reac_debug .and. i == reac_debug)then
                write(*,100) 'Forward rate constant is', lisol_work(1,1,1)
            endif
            ! forward rate of progress, mol/m3/s
            do j = 1,nspec
                lisol_work(:,:,:) = lisol_work(:,:,:) * lisol_props_species2(:,:,:,j)**(reac_coeffs(j,i))
            enddo
            ! debug
            if(if_reac_debug .and. i == reac_debug)then
                write(*,200) 'Forward rate of progress is', lisol_work(1,1,1)
            endif
            ! mass production rate, kg/m3/s
            do j = 1,nspec
                lisol_mass_prod_rate(:,:,:,j) = lisol_mass_prod_rate(:,:,:,j) &
                                                + devia_coeffs(j)*lisol_work(:,:,:)*lisol_molmas(j)
            enddo
            if(isRevers(i))then ! reversible reactions
                lisol_work2(:,:,:) = 0.d0
                do j = 1,nspec
                    lisol_work2(:,:,:) = lisol_work2(:,:,:) + devia_coeffs(j) * lisol_props_species(:,:,:,j)
                enddo
                rtemp = dble(sum(devia_coeffs))
                ! equilibrium constant
                lisol_work2(:,:,:) = dexp(-lisol_work2(:,:,:) / GasConstant / lisol_temperature(:,:,:) &
                                          + rtemp * lisol_heat_release(:,:,:)) ! temp array here
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,300) 'Equilibrium constant is', lisol_work2(1,1,1)
                endif
                ! back rate constant
                lisol_work2(:,:,:) = lisol_work3(:,:,:) / lisol_work2(:,:,:)
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,400) 'Back rate constant is', lisol_work2(1,1,1)
                endif
                ! back rate of progress, mol/m3/s
                do j = 1,nspec
                    lisol_work2(:,:,:) = lisol_work2(:,:,:) * lisol_props_species2(:,:,:,j)**(prod_coeffs(j,i))
                enddo
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,500) 'Back rate of progress is', lisol_work2(1,1,1)
                endif
                ! mass production rate, kg/m3/s
                do j = 1,nspec
                    lisol_mass_prod_rate(:,:,:,j) = lisol_mass_prod_rate(:,:,:,j) &
                                                    - devia_coeffs(j)*lisol_work2(:,:,:)*lisol_molmas(j)
                enddo
            endif
        ! 2. three-body reaction
        elseif(reac_type(i,1) .eq. THREE_BODY_RXN)then
            ! forward rate constant
            lisol_work(:,:,:) = rate_coeff(i,1) * lisol_temperature(:,:,:)**(rate_coeff(i,2)) &
                                * dexp(-rate_coeff(i,3) / GasConstant / lisol_temperature(:,:,:))
            lisol_work2(:,:,:) = lisol_work(:,:,:)
            ! debug
            if(if_reac_debug .and. i == reac_debug)then
                write(*,100) 'Forward rate constant is', lisol_work(1,1,1)
            endif
            ! forward rate of progress, mol/m3/s
            lisol_work3(:,:,:) = 0.d0
            do j = 1,nspec
                lisol_work(:,:,:) = lisol_work(:,:,:) * lisol_props_species2(:,:,:,j)**(reac_coeffs(j,i))
                lisol_work3(:,:,:) = lisol_work3(:,:,:) + collision_efficiency(i,j) * lisol_props_species2(:,:,:,j)
            enddo
            lisol_work(:,:,:) = lisol_work(:,:,:) * lisol_work3(:,:,:)
            ! debug
            if(if_reac_debug .and. i == reac_debug)then
                write(*,200) 'Forward rate of progress is', lisol_work(1,1,1)
            endif
            ! mass production rate, kg/m3/s
            do j = 1,nspec
                lisol_mass_prod_rate(:,:,:,j) = lisol_mass_prod_rate(:,:,:,j) &
                                                + devia_coeffs(j)*lisol_work(:,:,:)*lisol_molmas(j)
            enddo
            if(isRevers(i))then ! reversible reactions
                lisol_work(:,:,:) = 0.d0
                do j = 1,nspec
                    lisol_work(:,:,:) = lisol_work(:,:,:) + devia_coeffs(j) * lisol_props_species(:,:,:,j)
                enddo
                rtemp = dble(sum(devia_coeffs))
                ! equilibrium constant
                lisol_work(:,:,:) = dexp(-lisol_work(:,:,:) / GasConstant / lisol_temperature(:,:,:) &
                                          + rtemp * lisol_heat_release(:,:,:)) ! temp array here
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,300) 'Equilibrium constant is', lisol_work(1,1,1)
                endif
                ! back rate constant
                lisol_work2(:,:,:) = lisol_work2(:,:,:) / lisol_work(:,:,:)
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,400) 'Back rate constant is', lisol_work2(1,1,1)
                endif
                ! back rate of progress, mol/m3/s
                do j = 1,nspec
                    lisol_work2(:,:,:) = lisol_work2(:,:,:) * lisol_props_species2(:,:,:,j)**(prod_coeffs(j,i))
                enddo
                lisol_work2(:,:,:) = lisol_work2(:,:,:) * lisol_work3(:,:,:) ! collision efficiency
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,500) 'Back rate of progress is', lisol_work2(1,1,1)
                endif
                ! mass production rate, kg/m3/s
                do j = 1,nspec
                    lisol_mass_prod_rate(:,:,:,j) = lisol_mass_prod_rate(:,:,:,j) &
                                                    - devia_coeffs(j)*lisol_work2(:,:,:)*lisol_molmas(j)
                enddo
            endif
        ! 3. simple falloff reaction
        elseif((reac_type(i,1) .eq. FALLOFF_RXN) .and. (reac_type(i,2) .eq. SIMPLE_FALLOFF))then
            ! low pressure k0
            lisol_work(:,:,:) = rate_coeff(i,1) * lisol_temperature(:,:,:)**(rate_coeff(i,2)) &
                                * dexp(-rate_coeff(i,3) / GasConstant / lisol_temperature(:,:,:))
            ! high pressure, k_infty
            lisol_work2(:,:,:) = rate_coeff(i,4) * lisol_temperature(:,:,:)**(rate_coeff(i,5)) &
                                * dexp(-rate_coeff(i,6) / GasConstant / lisol_temperature(:,:,:))
            ! [M]
            lisol_work3(:,:,:) = 0.d0
            do j = 1,nspec
                lisol_work3(:,:,:) = lisol_work3(:,:,:) + collision_efficiency(i,j) * lisol_props_species2(:,:,:,j)
            enddo
            ! ! debug
            ! if(if_reac_debug .and. i == reac_debug)then
            !     write(*,*) '[M]', lisol_work3(1,1,1)
            ! endif
            ! reduced pressure Pr
            lisol_work4(:,:,:) = lisol_work3(:,:,:) * lisol_work(:,:,:) / lisol_work2(:,:,:)
            ! forward rate constant, Kf = K_infty * Pr / (1+Pr)
            lisol_work(:,:,:) = lisol_work2(:,:,:) * lisol_work4(:,:,:) / (1.d0 + lisol_work4(:,:,:))
            lisol_work3(:,:,:) = lisol_work(:,:,:)
            ! debug
            if(if_reac_debug .and. i == reac_debug)then
                write(*,100) 'Forward rate constant is', lisol_work(1,1,1)
            endif
            ! forward rate of progress, mol/m3/s
            do j = 1,nspec
                lisol_work(:,:,:) = lisol_work(:,:,:) * lisol_props_species2(:,:,:,j)**(reac_coeffs(j,i))
            enddo
            ! debug
            if(if_reac_debug .and. i == reac_debug)then
                write(*,200) 'Forward rate of progress is', lisol_work(1,1,1)
            endif
            ! mass production rate, kg/m3/s
            do j = 1,nspec
                lisol_mass_prod_rate(:,:,:,j) = lisol_mass_prod_rate(:,:,:,j) &
                                                + devia_coeffs(j)*lisol_work(:,:,:)*lisol_molmas(j)
            enddo
            if(isRevers(i))then ! reversible reactions
                lisol_work2(:,:,:) = 0.d0
                do j = 1,nspec
                    lisol_work2(:,:,:) = lisol_work2(:,:,:) + devia_coeffs(j) * lisol_props_species(:,:,:,j)
                enddo
                rtemp = dble(sum(devia_coeffs))
                ! equilibrium constant
                lisol_work2(:,:,:) = dexp(-lisol_work2(:,:,:) / GasConstant / lisol_temperature(:,:,:) &
                                          + rtemp * lisol_heat_release(:,:,:)) ! temp array here
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,300) 'Equilibrium constant is', lisol_work2(1,1,1)
                endif
                ! back rate constant
                lisol_work2(:,:,:) = lisol_work3(:,:,:) / lisol_work2(:,:,:)
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,400) 'Back rate constant is', lisol_work2(1,1,1)
                endif
                ! back rate of progress, mol/m3/s
                do j = 1,nspec
                    lisol_work2(:,:,:) = lisol_work2(:,:,:) * lisol_props_species2(:,:,:,j)**(prod_coeffs(j,i))
                enddo
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,500) 'Back rate of progress is', lisol_work2(1,1,1)
                endif
                ! mass production rate, kg/m3/s
                do j = 1,nspec
                    lisol_mass_prod_rate(:,:,:,j) = lisol_mass_prod_rate(:,:,:,j) &
                                                    - devia_coeffs(j)*lisol_work2(:,:,:)*lisol_molmas(j)
                enddo
            endif
        ! 4. TROE falloff reaction
        elseif((reac_type(i,1) .eq. FALLOFF_RXN) .and. (reac_type(i,2) .eq. TROE_FALLOFF))then
            ! low pressure k0
            lisol_work(:,:,:) = rate_coeff(i,1) * lisol_temperature(:,:,:)**(rate_coeff(i,2)) &
                                * dexp(-rate_coeff(i,3) / GasConstant / lisol_temperature(:,:,:))
            ! high pressure, k_infty
            lisol_work2(:,:,:) = rate_coeff(i,4) * lisol_temperature(:,:,:)**(rate_coeff(i,5)) &
                                * dexp(-rate_coeff(i,6) / GasConstant / lisol_temperature(:,:,:))
            lisol_work3(:,:,:) = 0.d0
            do j = 1,nspec
                lisol_work3(:,:,:) = lisol_work3(:,:,:) + collision_efficiency(i,j) * lisol_props_species2(:,:,:,j)
            enddo
            ! reduced pressure Pr
            lisol_work4(:,:,:) = lisol_work3(:,:,:) * lisol_work(:,:,:) / lisol_work2(:,:,:)
            ! K_infty * Pr / (1+Pr)
            lisol_work(:,:,:) = lisol_work2(:,:,:) * lisol_work4(:,:,:) / (1.d0 + lisol_work4(:,:,:))
            ! Fcent = (1-A)*exp(-T/T3) + A*exp(-T/T1) + exp(-T2/T), the last term is optional
            lisol_work2(:,:,:) = (1.d0 - rate_coeff(i,7)) * dexp(-lisol_temperature(:,:,:)/rate_coeff(i,8)) &
                                 + rate_coeff(i,7) * dexp(-lisol_temperature(:,:,:)/rate_coeff(i,9))
            if(rate_coeff(i,10) .ne. 0.d0)then
                lisol_work2(:,:,:) = lisol_work2(:,:,:) + dexp(-rate_coeff(i,10) / lisol_temperature(:,:,:)) 
            endif
            ! ! to avoid underflow warning
            ! flag = lisol_work2(:,:,:) > SmallNumber
            ! lisol_work2(:,:,:) = merge(lisol_work2(:,:,:), SmallNumber, flag)
            lisol_work2(:,:,:) = dlog10(lisol_work2(:,:,:))
            ! C = -0.4 - 0.67*log10(Fcent)
            lisol_work3(:,:,:) = -0.4d0 - 0.67d0 * lisol_work2(:,:,:)
            ! N = 0.75 - 1.27*log10(Fcent)
            lisol_work5(:,:,:) = 0.75d0 - 1.27d0 * lisol_work2(:,:,:)
            ! ! to avoid underflow warning
            ! flag = lisol_work4(:,:,:) > SmallNumber
            ! lisol_work4(:,:,:) = merge(lisol_work4(:,:,:), SmallNumber, flag)
            lisol_work4(:,:,:) = dlog10(lisol_work4(:,:,:))
            ! f1 = (log10(Pr) + C) / (N - 0.14*(log10(Pr) + C))
            lisol_work4(:,:,:) = (lisol_work4(:,:,:) + lisol_work3(:,:,:)) &
                                / (lisol_work5(:,:,:) - 0.14d0*(lisol_work4(:,:,:) + lisol_work3(:,:,:)))
            ! F = 10^(log10(Fcent) / (1+f1^2))
            lisol_work2(:,:,:) = 10.d0**(lisol_work2(:,:,:) / (1.d0 + lisol_work4(:,:,:)**(2.d0)))
            ! forward rate constant, Kf = Kf * F
            lisol_work(:,:,:) = lisol_work(:,:,:) * lisol_work2(:,:,:)
            lisol_work3(:,:,:) = lisol_work(:,:,:)
            ! debug
            if(if_reac_debug .and. i == reac_debug)then
                write(*,100) 'Forward rate constant is', lisol_work(1,1,1)
            endif
            ! forward rate of progress, mol/m3/s
            do j = 1,nspec
                lisol_work(:,:,:) = lisol_work(:,:,:) * lisol_props_species2(:,:,:,j)**(reac_coeffs(j,i))
            enddo
            ! debug
            if(if_reac_debug .and. i == reac_debug)then
                write(*,200) 'Forward rate of progress is', lisol_work(1,1,1)
            endif
            ! mass production rate, kg/m3/s
            do j = 1,nspec
                lisol_mass_prod_rate(:,:,:,j) = lisol_mass_prod_rate(:,:,:,j) &
                                                + devia_coeffs(j)*lisol_work(:,:,:)*lisol_molmas(j)
            enddo
            if(isRevers(i))then ! reversible reactions
                lisol_work2(:,:,:) = 0.d0
                do j = 1,nspec
                    lisol_work2(:,:,:) = lisol_work2(:,:,:) + devia_coeffs(j) * lisol_props_species(:,:,:,j)
                enddo
                rtemp = dble(sum(devia_coeffs))
                ! equilibrium constant
                lisol_work2(:,:,:) = dexp(-lisol_work2(:,:,:) / GasConstant / lisol_temperature(:,:,:) &
                                          + rtemp * lisol_heat_release(:,:,:)) ! temp array here
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,300) 'Equilibrium constant is', lisol_work2(1,1,1)
                endif
                ! back rate constant
                lisol_work2(:,:,:) = lisol_work3(:,:,:) / lisol_work2(:,:,:)
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,400) 'Back rate constant is', lisol_work2(1,1,1)
                endif
                ! back rate of progress, mol/m3/s
                do j = 1,nspec
                    lisol_work2(:,:,:) = lisol_work2(:,:,:) * lisol_props_species2(:,:,:,j)**(prod_coeffs(j,i))
                enddo
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,500) 'Back rate of progress is', lisol_work2(1,1,1)
                endif
                ! mass production rate, kg/m3/s
                do j = 1,nspec
                    lisol_mass_prod_rate(:,:,:,j) = lisol_mass_prod_rate(:,:,:,j) &
                                                    - devia_coeffs(j)*lisol_work2(:,:,:)*lisol_molmas(j)
                enddo
            endif
        ! 5. SRI falloff reaction
        elseif((reac_type(i,1) .eq. FALLOFF_RXN) .and. ((reac_type(i,2) .eq. SRI_FALLOFF_3) &
                                                   .or. (reac_type(i,2) .eq. SRI_FALLOFF_5)))then
            ! low pressure k0
            lisol_work(:,:,:) = rate_coeff(i,1) * lisol_temperature(:,:,:)**(rate_coeff(i,2)) &
                                * dexp(-rate_coeff(i,3) / GasConstant / lisol_temperature(:,:,:))
            ! high pressure, k_infty
            lisol_work2(:,:,:) = rate_coeff(i,4) * lisol_temperature(:,:,:)**(rate_coeff(i,5)) &
                                * dexp(-rate_coeff(i,6) / GasConstant / lisol_temperature(:,:,:))
            lisol_work3(:,:,:) = 0.d0
            do j = 1,nspec
                lisol_work3(:,:,:) = lisol_work3(:,:,:) + collision_efficiency(i,j) * lisol_props_species2(:,:,:,j)
            enddo
            ! reduced pressure Pr
            lisol_work4(:,:,:) = lisol_work3(:,:,:) * lisol_work(:,:,:) / lisol_work2(:,:,:)
            ! K_infty * Pr / (1+Pr)
            lisol_work(:,:,:) = lisol_work2(:,:,:) * lisol_work4(:,:,:) / (1.d0 + lisol_work4(:,:,:))
            ! ! to avoid underflow warning
            ! flag = lisol_work4(:,:,:) > SmallNumber
            ! lisol_work4(:,:,:) = merge(lisol_work4(:,:,:), SmallNumber, flag)
            ! F = d*[a*exp(-b/T) + exp(-T/c)]^[1/(1+(log10(Pr))^2)] * T^e
            lisol_work2(:,:,:) = rate_coeff(i,10) * (lisol_temperature(:,:,:)**(rate_coeff(i,11))) &
                        * (rate_coeff(i,7) * dexp(-rate_coeff(i,8)/lisol_temperature(:,:,:)) &
                         + dexp(-lisol_temperature(:,:,:)/rate_coeff(i,9)) )**(1.d0 / (1.d0 + (dlog10(lisol_work4(:,:,:)))**(2.d0)))
            ! forward rate constant, Kf = Kf * F
            lisol_work(:,:,:) = lisol_work(:,:,:) * lisol_work2(:,:,:)
            lisol_work3(:,:,:) = lisol_work(:,:,:)
            ! debug
            if(if_reac_debug .and. i == reac_debug)then
                write(*,100) 'Forward rate constant is', lisol_work(1,1,1)
            endif
            ! forward rate of progress, mol/m3/s
            do j = 1,nspec
                lisol_work(:,:,:) = lisol_work(:,:,:) * lisol_props_species2(:,:,:,j)**(reac_coeffs(j,i))
            enddo
            ! debug
            if(if_reac_debug .and. i == reac_debug)then
                write(*,200) 'Forward rate of progress is', lisol_work(1,1,1)
            endif
            ! mass production rate, kg/m3/s
            do j = 1,nspec
                lisol_mass_prod_rate(:,:,:,j) = lisol_mass_prod_rate(:,:,:,j) &
                                                + devia_coeffs(j)*lisol_work(:,:,:)*lisol_molmas(j)
            enddo
            if(isRevers(i))then ! reversible reactions
                lisol_work2(:,:,:) = 0.d0
                do j = 1,nspec
                    lisol_work2(:,:,:) = lisol_work2(:,:,:) + devia_coeffs(j) * lisol_props_species(:,:,:,j)
                enddo
                rtemp = dble(sum(devia_coeffs))
                ! equilibrium constant
                lisol_work2(:,:,:) = dexp(-lisol_work2(:,:,:) / GasConstant / lisol_temperature(:,:,:) &
                                          + rtemp * lisol_heat_release(:,:,:)) ! temp array here
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,300) 'Equilibrium constant is', lisol_work2(1,1,1)
                endif
                ! back rate constant
                lisol_work2(:,:,:) = lisol_work3(:,:,:) / lisol_work2(:,:,:)
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,400) 'Back rate constant is', lisol_work2(1,1,1)
                endif
                ! back rate of progress, mol/m3/s
                do j = 1,nspec
                    lisol_work2(:,:,:) = lisol_work2(:,:,:) * lisol_props_species2(:,:,:,j)**(prod_coeffs(j,i))
                enddo
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,500) 'Back rate of progress is', lisol_work2(1,1,1)
                endif
                ! mass production rate, kg/m3/s
                do j = 1,nspec
                    lisol_mass_prod_rate(:,:,:,j) = lisol_mass_prod_rate(:,:,:,j) &
                                                    - devia_coeffs(j)*lisol_work2(:,:,:)*lisol_molmas(j)
                enddo
                ! write(*,*) 'break point here'
            endif
        ! 6. PLOG reactions
        elseif(reac_type(i,1) .eq. PLOG_RXN)then
            do j = 1,Pnum_PLOG(i)-1
                where((rate_coeff_PLOG(i,j,1) <= lisol_pressure) .and. (lisol_pressure <= rate_coeff_PLOG(i,j+1,1)))
                    ! P1 -> k1(T)
                    lisol_work2 = rate_coeff_PLOG(i,j,2) * (lisol_temperature)**(rate_coeff_PLOG(i,j,3)) &
                                  * dexp(-rate_coeff_PLOG(i,j,4) / GasConstant / lisol_temperature)
                    lisol_work2 = dlog10(lisol_work2)
                    ! ! to avoid underflow warning
                    ! flag = lisol_work2 > SmallNumber
                    ! lisol_work2 = merge(lisol_work2, SmallNumber, flag)
                    ! P2 -> k2(T)
                    lisol_work3 = rate_coeff_PLOG(i,j+1,2) * (lisol_temperature)**(rate_coeff_PLOG(i,j+1,3)) &
                                  * dexp(-rate_coeff_PLOG(i,j+1,4) / GasConstant / lisol_temperature)
                    ! ! to avoid underflow warning
                    ! flag = lisol_work3 > SmallNumber
                    ! lisol_work3 = merge(lisol_work3, SmallNumber, flag)
                    ! forward rate constant
                    lisol_work = 10.d0**(lisol_work2 + (dlog10(lisol_work3) - lisol_work2) &
                                               * (dlog10(lisol_pressure) - dlog10(rate_coeff_PLOG(i,j,1))) &
                                               / (dlog10(rate_coeff_PLOG(i,j+1,1)) - dlog10(rate_coeff_PLOG(i,j,1))))
                end where
            enddo
            lisol_work3(:,:,:) = lisol_work(:,:,:)
            ! debug
            if(if_reac_debug .and. i == reac_debug)then
                write(*,100) 'Forward rate constant is', lisol_work(1,1,1)
            endif
            ! forward rate of progress, mol/m3/s
            do j = 1,nspec
                lisol_work(:,:,:) = lisol_work(:,:,:) * lisol_props_species2(:,:,:,j)**(reac_coeffs(j,i))
            enddo
            ! debug
            if(if_reac_debug .and. i == reac_debug)then
                write(*,200) 'Forward rate of progress is', lisol_work(1,1,1)
            endif
            ! mass production rate, kg/m3/s
            do j = 1,nspec
                lisol_mass_prod_rate(:,:,:,j) = lisol_mass_prod_rate(:,:,:,j) &
                                                + devia_coeffs(j)*lisol_work(:,:,:)*lisol_molmas(j)
            enddo
            if(isRevers(i))then ! reversible reactions
                lisol_work2(:,:,:) = 0.d0
                do j = 1,nspec
                    lisol_work2(:,:,:) = lisol_work2(:,:,:) + devia_coeffs(j) * lisol_props_species(:,:,:,j)
                enddo
                rtemp = dble(sum(devia_coeffs))
                ! equilibrium constant
                lisol_work2(:,:,:) = dexp(-lisol_work2(:,:,:) / GasConstant / lisol_temperature(:,:,:) &
                                          + rtemp * lisol_heat_release(:,:,:)) ! temp array here
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,300) 'Equilibrium constant is', lisol_work2(1,1,1)
                endif
                ! back rate constant
                lisol_work2(:,:,:) = lisol_work3(:,:,:) / lisol_work2(:,:,:)
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,400) 'Back rate constant is', lisol_work2(1,1,1)
                endif
                ! back rate of progress, mol/m3/s
                do j = 1,nspec
                    lisol_work2(:,:,:) = lisol_work2(:,:,:) * lisol_props_species2(:,:,:,j)**(prod_coeffs(j,i))
                enddo
                ! debug
                if(if_reac_debug .and. i == reac_debug)then
                    write(*,500) 'Back rate of progress is', lisol_work2(1,1,1)
                endif
                ! mass production rate, kg/m3/s
                do j = 1,nspec
                    lisol_mass_prod_rate(:,:,:,j) = lisol_mass_prod_rate(:,:,:,j) &
                                                    - devia_coeffs(j)*lisol_work2(:,:,:)*lisol_molmas(j)
                enddo
            endif
        else
            call lisol_error(1, &
            'lisol: error reac_type, not supported yet! Only for elementary, Three-body, Falloff, and PLOG reactions now.')
        endif
        !----
    enddo
    ! debug
    do i = 1,nspec
        if(spec_name(i)=='H2')then
            write(*,*) 'net production rate of '//trim(spec_name(i))//' is ', lisol_mass_prod_rate(1,1,1,i)
        endif
    enddo
    write(*,'(a)') '  Compute the temperature production rates'
    lisol_heat_prod_rate(:,:,:) = sum(lisol_enthalpy_spec(:,:,:,:) * lisol_mass_prod_rate(:,:,:,:), dim = 4)
    ! do i = 1,nspec
    !     lisol_heat_prod_rate(:,:,:) = lisol_heat_prod_rate(:,:,:) + lisol_enthalpy_spec(:,:,:,i) * lisol_mass_prod_rate(:,:,:,i)
    ! enddo
    ! heat release rate
    lisol_heat_release(:,:,:) = - lisol_heat_prod_rate(:,:,:)
    ! temperature source term
    lisol_heat_prod_rate(:,:,:) = - lisol_heat_prod_rate(:,:,:) / lisol_cp(:,:,:)
    ! debug
    write(*,'(a30,e24.16)') 'temperature production rate is', lisol_heat_prod_rate(1,1,1)
    write(*,'(a20,e24.16)') 'heat release rate is', lisol_heat_release(1,1,1)

100  format (a24,e24.16) ! Forward rate constant
200  format (a27,e24.16) ! Forward rate of progress
300  format (a23,e24.16) ! Equilibrium constant
400  format (a21,e24.16) ! Back rate constant
500  format (a24,e24.16) ! Back rate of progress
    return
end subroutine lisol_compute_source_terms
!=============================================
! This subroutine is used for forcing the sum of
! mass fraction equals to unity, and all fractions
! >=0 and <= 1
!=============================================
subroutine lisol_force_sum_massfrac_unity()
    use lisol_chemistry_mod
    implicit none
    integer(kind=4) :: l, flag_i
    if(y_forcing_type .eq. y_forcing_sum)then ! forcing the mass fraction by normalization
        write(*,'(a)') '  Forcing sum(Yi)=1 and all Yi >= 0 and <= 1; normalize all species'
        lisol_work(:,:,:) = sum(lisol_mass_frac(:,:,:,:), dim=4)
        do l = 1,nspec
            lisol_mass_frac(:,:,:,l) = lisol_mass_frac(:,:,:,l) / lisol_work(:,:,:)
        enddo
    elseif(y_forcing_type .eq. y_forcing_inert)then
        flag_i = 0
        do l = 1,nspec
            if(spec_name(l) .eq. trim(inert_gas))then
                flag_i = 1
                exit
            endif
        enddo
        if(flag_i .eq. 1)then
            write(*,'(a)') '  Forcing sum(Yi)=1 and all Yi >= 0 and <= 1; '//trim(inert_gas)//' as a garbage species'
            lisol_work(:,:,:) = 0.d0
            do l = 1,nspec
                if(spec_name(l) .ne. trim(inert_gas))then
                    lisol_work(:,:,:) = lisol_work(:,:,:) + lisol_mass_frac(:,:,:,l)
                endif
            enddo
            do l = 1,nspec
                if(spec_name(l) .eq. trim(inert_gas))then
                    lisol_mass_frac(:,:,:,l) = 1.d0 - lisol_work(:,:,:)
                endif
            enddo
        else
            call lisol_error(1, 'lisol_force_sum_massfrac_unity: error inert_gas '//trim(inert_gas)//', not found in the mechanism')
        endif
    endif

    return
end subroutine lisol_force_sum_massfrac_unity