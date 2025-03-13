        subroutine ch4_daycent(k)                         !outputs

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    estimate CH4 emission, calculating carbon substrate available for CH4 production, using carbon in active, slow, passive simulated from
!!    SWAT+ 
!!    ~ ~ ~ REFERENCE ~ ~ ~
!!    1. Hartman, M.D., et al., 2018, The Daily Century Ecosystem, Soil Organic Matter, Nutrient Cycling, Nitrogen Trace Gas, and 
!!       Methane Model User Manual, Scientific Basis, and Technical Documentation.
!!    2. Cheng, K., et al., 2013. Predicting methanogenesis from rice paddies using the DAYCENT ecosystem model. 
!!       Ecological Modelling 261–262, 19–31. https://doi.org/10.1016/j.ecolmodel.2013.04.003
!!    3. Huang, Y., et al., 2004. Modeling methane emission from rice paddies with various agricultural practices. 
!!       Journal of Geophysical Research: Atmospheres 109. https://doi.org/10.1029/2003JD004401
!!    4. Huang, YaO., et al., 1998. A semi-empirical model of methane emission from flooded rice paddy soils. 
!!       Global Change Biology 4, 247–268. https://doi.org/10.1046/j.1365-2486.1998.00129.x
!!    
!!    ~ ~ ~ NOTE ~ ~ ~
!!    1. In the methane simulation, the total soil organic carbon (SOC) from all SOM pools is aggregated and
!!       utilized as the carbon substrate for microbial processes driving methane production.
!!    2. Implement equations [1-eq26, 1-eq27] to to simulate soil redox potential (Eh) dynamics during flooding 
!!       and drainage events, using differential equations. 
!!    3. Verify whether root production refers to total root biomass or specifically root carbon production.
!!       Cross-check references [1-eq23] and [2-eq3] for clarification.
!!    4. Determine whether root production is expressed in kg biomass or kg C, and ensure values are 
!!       provided for each soil layer.
!!    5. The SWAT+ model includes a crop growth curve for plant development, seed dynamics, etc., 
!!       which can be used to estimate the maximum aboveground biomass at the end of the growing season.
!!    
    
    
        ! use hru_module, only : ihru
        use hru_module, only : ihru
        use soil_module
        use ghg_module
        use organic_mineral_mass_module
           
        implicit none
                
        integer :: j              !                     |number of hru
        integer,intent(in) :: k   !none                 |counter  
        
        ! local variable
        real :: eh_t                       ! mV             | Eh(t) represents the Eh value at time t, and t is the number 
                                           !                | of days after flooding began or since drainage occurred in the cycle. 
                                           !                | DEh and AEh are differential coefficients 
                                           !                | that were estimated as 0.16 and 0.23, respectively.
        real :: feh_t                      ! mV             | FEh, a reduction factor for soil redox potential (Eh) (mV), is estimated 
                                           !                | using the following equations from Huang et al. (1998) and Huang et al. (2004)
        
        real :: t_cont                     ! c'             | temperature threshold by Huang et al. (1998, 2004)
        real :: eh_cont                    ! mV             | Eh=-150mV for Eh<-150mV
        real :: f_t                        !                | temperature factor
        real :: subst_c                    ! kgC ha-1 day-1 | carbon subtrate
        real :: c_root                     ! kgC ha-1 day-1 | Root Carbon production (rhizodeposition)
        real :: si                         !                | SI is the soil texture index which is a function of the average sand content fraction         
        real :: sand_fr                    !                | (sand, 0.0 – 1.0) in the top 10 cm of soil  get from soil1
        real :: fp                         ! kg CH4-C ha-1 day-1 | FP is the fraction of CH4 emitted via rice plants.

        ! spark
        real :: ch4_prod                   ! kg CH4-C ha-1 day-1 | CH4 production rate
        real :: ch4ep                      ! kg CH4-C ha-1 day-1 | CH4 emissions through the rice plants
        real :: ch4ebl                     ! kg CH4-C ha-1 day-1 | Ebullition
        real :: ch4etot                    ! kg CH4-C ha-1 day-1 | Total CH4 emissions
        real :: max_abv_biom               ! kg biomass ha-1     | maximum aboveground biomass at the end of growing season,
                                           !                       theoretical maximum crop biomass where no more methane is emitted 
                                           !                       12600 (rice) value is provided in daycent.
        real :: root_mass                  ! kg biomass ha-1     | 
        ! spark        
        
         j = ihru
         
        ! hard user inputs
         eh_t = -100.0
         max_abv_biom = 200000. !
        ! root_mass = 1       

        ! soil texture index (SI) [1-eq22]
        sand_fr = 0.5 !spark need to get from soil properties (sand, silt, clay)
        si = 0.325 +2.25 * sand_fr         

        
        ! temperature function to get temp factor (q10, soil temperature) [1-eq24]
        if (soil(j)%phys(k)%tmp .gt. 30) then
          t_cont = 30.0
          f_t = dc_pars%q10 ** ((t_cont - 30.0) / 10.0)
        else
          f_t = dc_pars%q10 ** ((soil(j)%phys(k)%tmp - 30.0) / 10.0)
        end if
        
        ! c_substrate (kg ha-1 day-1)
        ! SOCs in Soil Organic Matter (SOM) pools are simulated using the "cbz_zhang2" subroutine. 
        ! In the methane simulation, the total soil organic carbon (SOC) from all SOM pools is aggregated and 
        ! utilized as the carbon substrate for microbial processes driving methane production.
        
        subst_c = soil1(j)%microb(k)%c + soil1(j)%meta(k)%c ! only active pools are considered
        ! subst_c = soil1(j)%microb(k)%c + soil1(j)%meta(k)%c + soil1(j)%hs(k)%c + soil1(j)%hp(k)%c    
          
        ! c_root: daily root carbon [1-eq23] 
        ! root_carbon should be considered by layers, but current variable is not.
        ! pl_mass(j)%root(k)%m represents root mass for each k layer
        ! pl_mass(j)%root_com%c should be the previous day's root production.
        
        ! Note:
        ! Verify whether we need root production or specifically root carbon production, check [2-eq3]
        c_root = dc_pars%fr_exudate * si * pl_mass(j)%root_com%c
        
        
        ! feh: estimate Eh factor [1-eq25] 
        ! we need to add an equation for getting eh_t for a flooding course [1-eq26] and
        ! drainage course [1-eq27].
        ! currently hard value -100 is used.
        ! check [1-eq26,27]
        if (eh_t >= -150) then
            feh_t = exp(-1.7*(150+eh_t)/150)
        else
            eh_cont = -150
            feh_t = exp(-1.7*(150+eh_cont)/150)
        end if
           
        
        ! ch4 production [1-eq20]
        ch4_prod = dc_pars%alpha1_ch4 * feh_t * (subst_c + (f_t * c_root)) ! layer specification required 
        
        !! ch4 emission
        ! emission through rice plant
        ! [1-eq29] uses a multiplier of 2.5 to convert carbon (C) to biomass. However, in this case, 
        ! the conversion is unnecessary because we directly provide aboveground biomass in kg biomass ha-1 instead.
        fp = 0.55 * (1 - ((pl_mass(j)%ab_gr_com%m / max_abv_biom))**0.25)
        ch4ep = fp * ch4_prod
        
        ! ebullition [1-eq30]
        if (pl_mass(j)%root_com%c .eq. 0) then
          ch4ebl = 0
        else
          ch4ebl = 0.7 * (ch4_prod - ch4ep- dc_pars%m0) * min(log(soil(j)%phys(k)%tmp), 1.0) * & 
                   (dc_pars%mrtblm / (pl_mass(j)%root_com%c * 2.5))
        end if
        
        ! ebullition can be estimated using [2-eq10] as well and root mass is required.
        ! pl_mass(j)%root(k)%m this is root mass
        ! following code produces error
        ! forrtl: severe (408): fort: (10): Subscript #1 of the array ROOT has value 4 
        ! which is greater than the upper bound of 3
        !if (pl_mass(j)%root(k)%m .eq. 0) then
        !  ch4ebl = 0
        !else
        !  ch4ebl = 0.7 * (ch4_prod - dc_pars%m0) * log(soil(j)%phys(k)%tmp / pl_mass(j)%root(k)%m) 
        !end if
        
        ch4etot = ch4ep + ch4ebl

        !write (4700,*) isd, time%day, time%yrc, pldb(iplt)%plantnm, hlt(isd)%alai, hlt(isd)%dm, yield
        print "(' hruid: ', i0, ' layer: ', i0, '  ch4_prod: ', es12.4, '  ch4ep: ', es12.4, ' ch4ebl: ', es12.4,'  ch4etot: ', es12.4)", &
            j, k, ch4_prod, ch4ep, ch4ebl, ch4etot
        return
        end subroutine ch4_daycent
          