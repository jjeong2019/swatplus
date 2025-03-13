        subroutine ch4_meres(k)                         !outputs

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    estimate CH4 emission, calculating carbon substrate available for CH4 production, using carbon in active, slow, passive simulated from
!!    SWAT+ 
!!    ~ ~ ~ REFERENCE ~ ~ ~
!!    1. Matthews et al., 2000. Using a Crop/Soil Simulation Model and GIS Techniques to Assess Methane Emissions from Rice Fields in 
!!       Asia. I. Model Development. Nutrient Cycling in Agroecosystems 58, 141–159. https://doi.org/10.1023/A:1009894619446
!!    
!!    ~ ~ ~ NOTE ~ ~ ~
!!    1. Except info SWAT+ can provide, theMERES model requires the following inputs: 
!!       O2, CH4, and AEA. These inputs are mandatory.
!!    2. Temperature and moisture effects are not considered in this model.
!!    3. The model requires root length density (m root m⁻³ soil) for each soil layer.

!!    ~ ~ ~ INPUT ~ ~ ~
        use hru_module, only : ihru
        use soil_module
        use ghg_module
        use organic_mineral_mass_module

        implicit none
                
        integer :: j              !                     |number of hru
        integer,intent(in) :: k   !none                 |counter  
        
        ! local variable (kg/ha/d=(mol C m⁻³ d⁻¹)×120.1×D)
        ! spark
        real :: subst_c                     ! kgC ha-1 day-1 | carbon subtrate
        real :: aea_ox                      ! kg CH4-C ha-1 day-1 | CH4 production rate
        real :: aea_ox_c                    ! mol C m-3           | the critical concentration of the oxidized alternative electron acceptor pool
        real :: pch4prod                    ! kg CH4-C ha-1 day-1 | potential CH4 production        
        real :: ach4prod                    ! kg CH4-C ha-1 day-1 | actual CH4 production        
        real :: sdepth                      ! meters              | depth of soil layer
        real :: chg_aea_ox                  ! meters              | the critical concentration of the oxidized alternative electron acceptor pool 
        real :: o2conc                      ! mol m-3             | the concentration of O2
        real :: ch4conc                     ! mol m-3             | the concentration of CH4
        real :: ch4cons                     ! kg CH4-C ha-1 day-1 | the concentration of CH4
        real :: o2cons                      ! mol m-3             | the concentration of CH4
        real :: o2flux                      ! mol m-3             | the concentration of CH4
        real :: ch4flux                      ! mol m-3             | the concentration of CH4
        !

        ! spark        
        j = ihru
        
        ! hard user inputs
        aea_ox = 100. ! mol C m-3
        o2conc = 100. ! mol m-3
        ch4conc = 100. ! mol m-3

        ! The effect of alternative electron acceptors on CH4 production (AEA)
        ! CH2O + AEAox (oxidized form) -> CO2 + AEAred (reduced) + 2H+
        sdepth = soil(j)%phys(k)%d *0.001 ! meters
        
        ! carbon substrate from SWAT+
        ! kgC ha-1 day-1 = kgC ha-1 day-1 + kgC ha-1 day-1
        subst_c = soil1(j)%microb(k)%c + soil1(j)%meta(k)%c ! only active pools are considered

        ! potential methane production [1-eq7]
        ! aea_ox_c critical, can be parameterized 
        aea_ox_c = 24.0 * 120.1 * sdepth ! mol C m-3 to kg/ha/d
        if (aea_ox .gt. aea_ox_c) then
          pch4prod = 0.0
        else if ((aea_ox .gt. 0.0) .and. (aea_ox .lt. aea_ox_c)) then
          pch4prod = min(0.2* (1-(aea_ox/aea_ox_c)), subst_c)
        else if  (aea_ox .eq. 0.0) then
          pch4prod = subst_c
        end if
        
        ! the rage of change of the oxidized aea pool [1-eq8, kg/ha/d]
        chg_aea_ox = subst_c - 2.0 * pch4prod

        ! reoxidation of the aea pool in midseason drainage
        ! Actual CH4 production (CH4) in a given soil layer [1-eq8, kg/ha/d]
        ! kg/ha/d = kg/ha/d / (1.0 + (mol/m3 * mol m-3))
        ach4prod = pch4prod / (1.0 + (meres_pars%eta * o2conc)) ! meres_pars%eta * o2conc will be a ratio

        ! The rate of CH4 consumption (QCH4, kg CH4-C ha-1 day-1) by 
        ! the methanotrophic bacteria (see 1-eq13) in 
        ! a soil layer is given by the Michaelis-Menten equation
        ! kg/ha/d = kg/ha/d * (mol/m3 / (mol/m3 + mol/m3)) * (mol/m3 / (mol/m3 + mol/m3))
        ch4cons = (pch4prod *                                          &
            (ch4conc/(meres_pars%k1 + ch4conc)) *               &
            (o2conc/(meres_pars%k2 + o2conc)))

            
        ! Oxygen consumption rate (QO2, kg CH4-C ha-1 day-1) [1-eq14]
        ! kg/ha/d = kg/ha/d * (mol/m3 / (mol/m3 + mol/m3)) 
        o2cons = 2*ch4cons + 2*pch4prod * (o2conc / (meres_pars%k3 + o2conc))


        ! we need to calculate root length density (m root m-3 soil), rld 
        ! the flux (F, mol m-3 s-1) for each substance (O2 or CH4) is given by [1-eq15]
        ! convert seconds to day
        ! A positive value for F represents flux from atmosphere to soil, and a negative value vice versa.
        ! mol/m3/d = m/m * (m root m-3 soil) * m2/d * (mol/m3 - mol/m3)

        o2flux = (                                            &
        meres_pars%sc_root * (0.1 * 10e+04) *                 &
        meres_pars%o2_dif * (meres_pars%o2_ya0 - o2conc)      &
        )        

        ch4flux = (                                            &
            meres_pars%sc_root * (0.1 * 10e+04) *              &
            meres_pars%ch4_dif * (meres_pars%ch4_ya0 - ch4conc) &
        )        

        print "(' hruid: ', i0, ' layer: ', i0, '  depth: ', es12.4, '  chpp: ', es12.4 '  ch4flux: ', es12.4)", &
            j, k, sdepth, pch4prod, ch4flux

        return
        end subroutine ch4_meres
          