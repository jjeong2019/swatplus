       subroutine ch4_dndc(k)                         !outputs

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    estimate CH4 emission, calculating carbon substrate available for CH4 production, using carbon in active, slow, passive simulated from
!!    SWAT+ 
!!    ~ ~ ~ REFERENCE ~ ~ ~
!!    1. Li, C. S. (2000), Modeling trace gas emissions from agricultural ecosystems, Nutrient Cycling Agroecosystem, 
!!       58 (1–3), 259–276, doi:10.1023/ A:1009859006242
!!    2. DNDC_Scientific_Basis_and_Processes, 2017. . 
!!       Institute for the Study of Earth, Oceans, and Space University of New Hampshire, Durham, NH 03824, USA.

!!    
!!    ~ ~ ~ NOTE ~ ~ ~
!!    1. no root exudation is considered
!!    2. plant growth index (PGI) is required, days after planting (DAP) is used to calculate PGI
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
       !spark
       real :: poro                       !                     | soil porosity by layer
       real :: temp                       !                     | soil temperature by layer     
       real :: subst_c                    ! kgC ha-1 day-1 | carbon subtrate
       
       real :: ch4p                       ! kg CH4-C ha-1 day-1 | CH4 production rate
       real :: f_t                        !                     | temperature factor for 1-eq3.1
       real :: pgi                        !                     | plant growth index        
       real :: aere                       !                     | plant aerenchyma
       real :: ch4aere                    ! kg CH4-C ha-1 day-1 | Methane flux through plant aerenchyma
       
       real :: ch4eb                      ! kg CH4-C ha-1 day-1 | Methane flux through ebullition
       real :: ch4difu                    ! kg CH4-C ha-1 day-1 | Methane diffusion   
       
       real :: ch4o                       ! kg O-C ha-1 day-1   | Methane oxidation rate
       real :: f_t2                       !                     | temperature factor for 1-eq3.4        
       
       ! Institute for the Study of Earth, Oceans, and Space 2017 
       real :: oxidant                    ! mole/liter          | concentration of the dominant oxidant in the system        
       real :: reductant                  ! mole/liter          | concentration of the dominant reductant in the system           
       real :: eh                         ! V                   | redox potential of an oxidation-reduction system
       real :: e0                         ! V                   | standard electromotive force   
       integer :: nelectron               ! count               | the transferred electron number
       real :: eh2                        ! V                   | If we set temperature and pH to be constants, 
                                                               ! the Nernst equation can be simplified as (Li et al. 1994)       
                                                               ! (pH 7.0, temperature 25°C)          
       real :: o2                         ! mole/liter          | oxygen concentration         
       
       real :: ch4etot                    ! kg CH4-C ha-1 day-1 | Total CH4 emissions


       ! spark        
       
       j = ihru
       

       ! hard user inputs
       pgi = 1./90.
       oxidant = 1.
       reductant = 100.
       e0 = 1.
       nelectron = 100
       o2 = 100.
       
       
       
       
       
       !eh_t = -100.0
       !max_abv_biom = 200000. !
       
       ! root_mass = 1       

       poro = soil(j)%phys(k)%por
       temp = soil(j)%phys(k)%tmp

       
       ! temperature function to get temp factor [1-eq3.1]
       f_t = dndc_pars%b * exp(0.2424 * temp)
       
       ! temperature function to get temp factor [1-eq3.4]
       f_t2 = (-0.1687 * (0.1 * temp)**3) + (1.167 * (0.1 * temp)**2 ) - &
              (2.0303 * (0.1 * temp) + 1.042)
       
       ! c_substrate (kg ha-1 day-1)
       ! SOCs in Soil Organic Matter (SOM) pools are simulated using the "cbz_zhang2" subroutine. 
       ! In the methane simulation, the total soil organic carbon (SOC) from all SOM pools is aggregated and 
       ! utilized as the carbon substrate for microbial processes driving methane production.
       
       subst_c = soil1(j)%microb(k)%c + soil1(j)%meta(k)%c ! only active pools are considered 
       ! subst_c = soil1(j)%microb(k)%c + soil1(j)%meta(k)%c + soil1(j)%hs(k)%c + soil1(j)%hp(k)%c    
              
       
       
       
       ! eh: redox potential  
       eh = e0 + ((dndc_pars%Rgas * (temp+273.)) / (nelectron * dndc_pars%Faraday)) * log(oxidant / reductant)
       
       ! eh2: (pH 7.0, temperature 25°C)
       eh2 = 0.82 + 0.015 * log(o2)
       
       
       
       ! ch4 production [1-eq3.1]
       ch4p = dndc_pars%a * subst_c * f_t  ! AC is subst_c
       
       !! ch4 emission
       ! emission through plant aerenchyma
       pgi = 30. / 90.
       aere = (-0.0009 * pgi**5) + (0.0047 * pgi**4) + (-0.883 * pgi**3) + (1.9863 * pgi**2)  &
              + (-0.3795 * pgi +0.0251)
       ch4aere = 0.5 * ch4p * aere

       
       ! ebullition [1-eq3.4]
       ch4eb = 0.025 * ch4p* poro * f_t2 * (1 - aere)
       
       ! diffusion [1-eq3.5]
       ch4difu = 0.01 * ch4p - ch4p * soil(j)%phys(k)%tmp * poro
       

       ch4etot = ch4aere + ch4eb + ch4difu
       !
       !write (4700,*) isd, time%day, time%yrc, pldb(iplt)%plantnm, hlt(isd)%alai, hlt(isd)%dm, yield
       print "(' hruid: ', i0, ' layer: ', i0, '  ch4_prod: ', es12.4, '  ch4ep: ', es12.4, ' ch4ebl: ', es12.4,'  ch4etot: ', es12.4)", &
       j, k, ch4p, ch4aere, ch4eb, ch4etot
       return
       end subroutine ch4_dndc
       