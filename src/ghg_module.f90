      module ghg_module 
      
      implicit none

      !! The diffusion of CH4 across the air–water interface is a relatively minor pathway for transporting CH4 to the atmosphere 
      !! (Shearer and Khalil, 1993; van Bodegom et al., 2001) and is not addressed in the DAYCENT model. 
      !! However, we plan to include the diffusion function in the development of the GHG module in SWAT+.

      ! DayCent
      type dc_parms
        ! decomposition per day Dangal et al., 2022. 
        ! Improving Soil Carbon Estimates by Linking Conceptual Pools Against Measurable Carbon Fractions in 
        ! the DAYCENT Model Version 4.5. Journal of Advances in Modeling Earth Systems 14, e2021MS002622. 
        ! https://doi.org/10.1029/2021MS002622
          
        real :: q10 = 3                     !           | temperature coefficient representing the change of a biological or 
                                            !           | chemical system as a consequence of increasing the temperature by 10°C, 
                                            !           | and was assumed to be a value of 3.0 (Huang et al. 1998). Defaults to 3.0.
        real :: fr_exudate = 0.45           !           | fraction of root carbon production contributing to rhizodeposition (fracToExudates)
                                            !           | range 0.35-0.6 beta2 FREXUD 0.0 - 1.0
        real :: alpha1_ch4 = 0.5            !           | conversion factor of carbohydrate decomposition to CO2 and CH4 (dimensionless) set to 0.5 [1]
        real :: b2_ch4 = 0.45               !           | b2 is the fraction of root carbon production contributing to rhizodeposition, 
                                            !           | which was set to an intermediate value of 0.45 from a range of 0.35–0.6 described 
                                            !           | in Cao et al. (1995). This range was based on previous studies 
                                            !           |  (Lambers, 1987; Merckx et al., 1985; Keith et al., 1986; Martin and Kemp, 1986; 
                                            !           | Buyanovsky and Wagner, 1987).
        real :: m0  = 0.015                 ! kgC ha−2 d−1    | Mo is set to 0.0015 gC m−2 d−1 (Huang et al., 2004),  but is set to 0.002 in DayCent           
        real :: mrtblm  = 10.0              ! kg biomass ha−2 | Root biomass that when exceeded start to reduce methane bubble formation
      end type dc_parms
      type (dc_parms) :: dc_pars
      
      ! DNDC
      type dndc_parms
        real :: a = 0.47                  !                 | This constant has been determined through empirical studies and model calibration.
        real :: b = 0.0009                !                 | This constant has been determined through empirical studies and model calibration.
        real :: Rgas = 8.314              ! 8.314 J/mol/k   | gas constant
        real :: Faraday = 96485           ! c/mol/          | Faraday constant       
      end type dndc_parms
      type (dndc_parms) :: dndc_pars
   
      ! MERES
      type meres_parms
        real :: aea_ox_c                   ! mol C m-3          | the critical concentration of the oxidized alternative electron acceptor pool
        real :: eta = 400                  ! m3 mol-1           | representing the sensitivity of methanogenesis to the concentration of O2
                                                                ! A value of 400 m3 mol-1 was used for η (Arah & Stephen, 1998)., defaults to 400
        ! k1, k2 and k3 are Michaelis-Menten constants for a dual-substrate reaction.
        real :: k1 = 0.33                  ! mol m-3            | k1: 
        real :: k2 = 0.44                  ! mol m-3            |             
        real :: k3 = 0.22                  ! mol m-3            |            
        real :: sc_root = 3.0e-04          ! m air (m root)-1   | the specific conductivity of the root system., defaults to 3.0e-04           
        real :: ch4_ya0 = 7.76             ! mol m-3            | ch4_ya0 is the concentration of the ch4 substance in the atmosphere           
        real :: o2_ya0 = 7.5e-05           ! mol m-3            | o2_ya0 is the concentration of the O2 substance in the atmosphere           
        real :: o2_dif = 1.74528           ! m2/d               | Diffusion constants (Da) of O2 (2.02 × 10-5 m2 s-1)           
        real :: ch4_dif = 0.91584          ! m2/d               | Diffusion constants (Da) of ch4 (1.06 × 10-5 m2 s-1)           
      end type meres_parms
      type (meres_parms) :: meres_pars      

      end module ghg_module