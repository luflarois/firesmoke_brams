 MODULE mod_chem_spack_kinetic
  
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: kinetic ! subroutine
 CONTAINS
   SUBROUTINE kinetic(jppj,Jphoto,rk,temp,xlw,Press,cosz,att,ijkbeg,ijkend,maxblock_size,nr)
 
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the kinetic rates for the gas-phase.
!     This routine is automatically generated by SPACK.
!     Mechanism: ../Mechanism/RELACS 
!     Species: ../Mechanism/ciRLCS 
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     TEMP: temperature ([K]).
!     XLW: water massic fraction.
!     PRESS: pressure ([Pa]).
!     ATT: attenuation variable.
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     RK: kinetic rates.
!
!------------------------------------------------------------------------
!
!     -- REMARKS
!
!------------------------------------------------------------------------
!
!     -- MODIFICATIONS
!
!------------------------------------------------------------------------
!
!     -- AUTHOR(S)
!
!     SPACK.
!
!------------------------------------------------------------------------
 
      IMPLICIT NONE
 
 
 
      INTEGER,INTENT(IN) :: jppj,ijkbeg,ijkend,maxblock_size,nr
      DOUBLE PRECISION,INTENT(IN) ::  xlw(maxblock_size),att(maxblock_size),cosz(maxblock_size)
      DOUBLE PRECISION,INTENT(IN) ::  temp(maxblock_size),Press(maxblock_size)
      DOUBLE PRECISION,INTENT(OUT) ::  rk(maxblock_size,nr)
      DOUBLE PRECISION,DIMENSION(maxblock_size) ::  Effko,Rapk,facteur
      DOUBLE PRECISION,DIMENSION(maxblock_size) :: YlH2O,SumM,azi
      DOUBLE PRECISION,INTENT(IN) ::  Jphoto(maxblock_size,jppj)
      INTEGER :: ijk
 
!     Compute third body.
!     Conversion = Avogadro*1d-6/Perfect gas constant.
!     PRESS in Pascal, SUMM in molecules/cm3, TEMP in Kelvin
 
      !DO ijk=ijkbeg,ijkend
        SumM(ijkbeg:ijkEnd) = Press(ijkbeg:ijkend) * 7.243D16 / temp(ijkbeg:ijkend)
      !END DO
 
!     Number of water molecules computed from the massic fraction
!     (absolute humidity)
 
     ! DO ijk=ijkbeg,ijkend
         YlH2O(ijkbeg:ijkend) = 29.d0*SumM(ijkbeg:ijkend)*xlw(ijkbeg:ijkend)/(18.d0+11.d0*xlw(ijkbeg:ijkend))
     ! END DO
 
!     For the zenithal angle at tropics
 
 
 
     !  DO ijk=ijkbeg,ijkend
       rk(ijkbeg:ijkend,  1) = Jphoto(ijkbeg:ijkend,  1)
       rk(ijkbeg:ijkend,  2) = Jphoto(ijkbeg:ijkend,  2)
       rk(ijkbeg:ijkend,  3) = Jphoto(ijkbeg:ijkend,  3)
       rk(ijkbeg:ijkend,  4) = Jphoto(ijkbeg:ijkend,  4)
       rk(ijkbeg:ijkend,  5) = Jphoto(ijkbeg:ijkend,  5)
       rk(ijkbeg:ijkend,  6) = Jphoto(ijkbeg:ijkend,  6)
       rk(ijkbeg:ijkend,  7) = Jphoto(ijkbeg:ijkend,  7)
       rk(ijkbeg:ijkend,  8) = Jphoto(ijkbeg:ijkend,  8)
       rk(ijkbeg:ijkend,  9) = Jphoto(ijkbeg:ijkend,  9)
       rk(ijkbeg:ijkend, 10) = Jphoto(ijkbeg:ijkend, 10)
       rk(ijkbeg:ijkend, 11) = Jphoto(ijkbeg:ijkend, 11)
       rk(ijkbeg:ijkend, 12) = Jphoto(ijkbeg:ijkend, 12)
       rk(ijkbeg:ijkend, 13) = Jphoto(ijkbeg:ijkend, 13)
       rk(ijkbeg:ijkend, 14) = Jphoto(ijkbeg:ijkend, 14)
       rk(ijkbeg:ijkend, 15) = Jphoto(ijkbeg:ijkend, 15)
       rk(ijkbeg:ijkend, 16) = Jphoto(ijkbeg:ijkend, 16)
       rk(ijkbeg:ijkend, 17) = Jphoto(ijkbeg:ijkend, 17)
      rk(ijkbeg:ijkend, 18) = SumM(ijkbeg:ijkend) * 6.0d-34 * (temp(ijkbeg:ijkend)/3.d2) ** (-2.3d0)
      rk(ijkbeg:ijkend, 18) = rk(ijkbeg:ijkend, 18) * SumM(ijkbeg:ijkend) * 0.2d0
      rk(ijkbeg:ijkend, 19) =  DEXP(-0.2555157957424871D+02   &
         - (  0.2060000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 20) =  DEXP(-0.2474064935803238D+02   &
         - ( -0.1100000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 20) = rk(ijkbeg:ijkend, 20) * SumM(ijkbeg:ijkend) * 0.8d0
      rk(ijkbeg:ijkend, 21) =  DEXP(-0.2416528521312882D+02   &
         - ( -0.7000000000000000D+02 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 21) = rk(ijkbeg:ijkend, 21) * SumM(ijkbeg:ijkend) * 0.2d0
      rk(ijkbeg:ijkend, 22) =  0.2200000000000000D-09
      rk(ijkbeg:ijkend, 22) = rk(ijkbeg:ijkend, 22) * YlH2O(ijkbeg:ijkend)
      rk(ijkbeg:ijkend, 23) =  DEXP(-0.2716101748668281D+02   &
         - (  0.9400000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 24) =  DEXP(-0.3214088112211231D+02   &
         - (  0.5000000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 25) =  DEXP(-0.2375982010502066D+02   &
         - ( -0.2500000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 26) =  DEXP(-0.2656631037893612D+02   &
         - (  0.1600000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 27) = 2.2d-13 * dexp(620.0d0 / temp(ijkbeg:ijkend))   &
                    + 1.9d-33* SumM(ijkbeg:ijkend) * dexp(980.0d0 / temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 28) = 3.08d-34 * dexp(2820.0d0 / temp(ijkbeg:ijkend)) +    &
                    2.66d-54 * SumM(ijkbeg:ijkend) * dexp(3180.0d0 / temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 28) = rk(ijkbeg:ijkend, 28) * YlH2O(ijkbeg:ijkend)
      Effko(ijkbeg:ijkend) =  0.9000000000000000D-31* (temp(ijkbeg:ijkend) / 3.d2)   &
                   **(- ( 0.1500000000000000D+01))
      Rapk(ijkbeg:ijkend) =  0.3000000000000000D-10* (temp(ijkbeg:ijkend) / 3.d2)   &
                    **(- ( 0.0000000000000000D+00))
      rk(ijkbeg:ijkend, 29) = (Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) /    &
                    ( 1.0d0 + Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend)))**2))
      rk(ijkbeg:ijkend, 30) =  DEXP(-0.2575921893902696D+02   &
         - ( -0.1200000000000000D+03 )/temp(ijkbeg:ijkend))
      Effko(ijkbeg:ijkend) =  0.9000000000000000D-31* (temp(ijkbeg:ijkend) / 3.d2)   &
                   **(- ( 0.2000000000000000D+01))
      Rapk(ijkbeg:ijkend) =  0.2200000000000000D-10* (temp(ijkbeg:ijkend) / 3.d2)   &
                    **(- ( 0.0000000000000000D+00))
      rk(ijkbeg:ijkend, 31) = (Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) /    &
                    ( 1.0d0 + Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend)))**2))
      Effko(ijkbeg:ijkend) =  0.7000000000000000D-30* (temp(ijkbeg:ijkend) / 3.d2)   &
                   **(- ( 0.2600000000000000D+01))
      Rapk(ijkbeg:ijkend) =  0.1500000000000000D-10* (temp(ijkbeg:ijkend) / 3.d2)   &
                    **(- ( 0.5000000000000000D+00))
      rk(ijkbeg:ijkend, 32) = (Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) /    &
                    ( 1.0d0 + Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend)))**2))
      Effko(ijkbeg:ijkend) =  0.2600000000000000D-29* (temp(ijkbeg:ijkend) / 3.d2)   &
                   **(- ( 0.3200000000000000D+01))
      Rapk(ijkbeg:ijkend) =  0.2400000000000000D-10* (temp(ijkbeg:ijkend) / 3.d2)   &
                    **(- ( 0.1300000000000000D+01))
      rk(ijkbeg:ijkend, 33) = (Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) /    &
                    ( 1.0d0 + Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend)))**2))
      rk(ijkbeg:ijkend, 34) =  0.2200000000000000D-10
      rk(ijkbeg:ijkend, 35) =  DEXP(-0.2632268829627837D+02   &
         - ( -0.2500000000000000D+03 )/temp(ijkbeg:ijkend))
      Effko(ijkbeg:ijkend) =  0.1800000000000000D-30* (temp(ijkbeg:ijkend) / 3.d2)   &
                   **(- ( 0.3200000000000000D+01))
      Rapk(ijkbeg:ijkend) =  0.4700000000000000D-11* (temp(ijkbeg:ijkend) / 3.d2)   &
                    **(- ( 0.1400000000000000D+01))
      rk(ijkbeg:ijkend, 36) = (Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) /    &
                    ( 1.0d0 + Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend)))**2))
      rk(ijkbeg:ijkend, 37) =  DEXP( 0.6142746008608852D+02   &
         - (  0.1090000000000000D+05 )/temp(ijkbeg:ijkend))
      Effko(ijkbeg:ijkend) =  0.1800000000000000D-30* (temp(ijkbeg:ijkend) / 3.d2)   &
                   **(- ( 0.3200000000000000D+01))
      Rapk(ijkbeg:ijkend) =  0.4700000000000000D-11* (temp(ijkbeg:ijkend) / 3.d2)   &
                    **(- ( 0.1400000000000000D+01))
      facteur(ijkbeg:ijkend) = (Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / ( 1.0d0 + Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) /    &
                 Rapk(ijkbeg:ijkend))) * 0.6d0 ** (1.0d0 / (1.0d0 +    &
                  (LOG10(Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend)))**2))
      rk(ijkbeg:ijkend, 37) = facteur(ijkbeg:ijkend) * rk(ijkbeg:ijkend, 37)
      rk(ijkbeg:ijkend, 38) =  0.3500000000000000D-11
      rk(ijkbeg:ijkend, 39) =  DEXP(-0.2474064935803238D+02   &
         - (  0.3900000000000000D+03 )/temp(ijkbeg:ijkend))
      Effko(ijkbeg:ijkend) = 7.2d-15 * dexp(785.0d0 / temp(ijkbeg:ijkend))
      Rapk(ijkbeg:ijkend) = 4.1d-16 * dexp(1440.0d0 / temp(ijkbeg:ijkend))
      facteur(ijkbeg:ijkend) =1.9d-33 * dexp(725.0d0 / temp(ijkbeg:ijkend)) * SumM(ijkbeg:ijkend)
      rk(ijkbeg:ijkend, 40) = Effko(ijkbeg:ijkend) + facteur(ijkbeg:ijkend)/(1.0d0 + facteur(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 41) =  DEXP(-0.2736865685146106D+02   &
         - ( -0.3800000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 42) =  DEXP(-0.2693787393536860D+02   &
         - (  0.1400000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 43) =  DEXP(-0.2975128465212864D+02   &
         - (  0.2450000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 44) =  0.3300000000000000D-38
      rk(ijkbeg:ijkend, 44) = rk(ijkbeg:ijkend, 44) * SumM(ijkbeg:ijkend) * 0.2d0
      rk(ijkbeg:ijkend, 45) =  DEXP(-0.2492297091482634D+02   &
         - ( -0.1700000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 46) =  DEXP(-0.3073211390514037D+02   &
         - (  0.1260000000000000D+04 )/temp(ijkbeg:ijkend))
      Effko(ijkbeg:ijkend) =  0.2200000000000000D-29* (temp(ijkbeg:ijkend) / 3.d2)   &
                   **(- ( 0.3900000000000000D+01))
      Rapk(ijkbeg:ijkend) =  0.1500000000000000D-11* (temp(ijkbeg:ijkend) / 3.d2)   &
                    **(- ( 0.7000000000000000D+00))
      rk(ijkbeg:ijkend, 47) = (Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) /    &
                    ( 1.0d0 + Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend)))**2))
      rk(ijkbeg:ijkend, 48) =  DEXP( 0.6117554523749536D+02   &
         - (  0.1100000000000000D+05 )/temp(ijkbeg:ijkend))
      Effko(ijkbeg:ijkend) =  0.2200000000000000D-29* (temp(ijkbeg:ijkend) / 3.d2)   &
                   **(- ( 0.3900000000000000D+01))
      Rapk(ijkbeg:ijkend) =  0.1500000000000000D-11* (temp(ijkbeg:ijkend) / 3.d2)   &
                    **(- ( 0.7000000000000000D+00))
      facteur(ijkbeg:ijkend) = (Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / ( 1.0d0 + Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) /    &
                 Rapk(ijkbeg:ijkend))) * 0.6d0 ** (1.0d0 / (1.0d0 +    &
                  (LOG10(Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend)))**2))
      rk(ijkbeg:ijkend, 48) = facteur(ijkbeg:ijkend) * rk(ijkbeg:ijkend, 48)
      rk(ijkbeg:ijkend, 49) =  DEXP(-0.2779354004542632D+02   &
         - (  0.2450000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 50) =  DEXP(-0.2592627302369012D+02   &
         - (  0.2000000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 50) = rk(ijkbeg:ijkend, 50) * SumM(ijkbeg:ijkend) * 5.8d-7
      Effko(ijkbeg:ijkend) =  0.3000000000000000D-30* (temp(ijkbeg:ijkend) / 3.d2)   &
                   **(- ( 0.3300000000000000D+01))
      Rapk(ijkbeg:ijkend) =  0.1500000000000000D-11* (temp(ijkbeg:ijkend) / 3.d2)   &
                    **(- ( 0.0000000000000000D+00))
      rk(ijkbeg:ijkend, 51) = (Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) /    &
                    ( 1.0d0 + Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend)))**2))
      rk(ijkbeg:ijkend, 52) = 1.5d-13 * (1.0d0 + 2.439d-20 * SumM(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 53) =  0.6000000000000000D-10
      rk(ijkbeg:ijkend, 54) =  0.0000000000000000D+00
      rk(ijkbeg:ijkend, 55) =  DEXP(-0.3943966082504782D+02   &
        + ( 0.2000000000000000D+01 * LOG(temp(ijkbeg:ijkend)))   &
        -  0.1361000000000000D+04/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 56) =  DEXP(-0.3873183693007194D+02   &
        + ( 0.2000000000000000D+01 * LOG(temp(ijkbeg:ijkend)))   &
        -  0.4920000000000000D+03/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 57) =   0.3760000000000000D-11   &
         * DEXP(-(  0.2600000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.1700000000000000D-11   &
         * DEXP(-(  0.1550000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.1210000000000000D-11   &
         * DEXP(-(  0.1250000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 58) =   0.1780000000000000D-11   &
         * DEXP(-( -0.4380000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.6070000000000000D-12   &
         * DEXP(-( -0.5000000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.0000000000000000D+00   &
         * DEXP(-( -0.4480000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 59) =  0.2540000000000000D-10   &
         * DEXP(-( -0.4100000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.0000000000000000D+00   &
         * DEXP(-( -0.4400000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 60) =  0.3310000000000000D-11   &
        * DEXP(-( -0.3550000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.3450000000000000D-12
      rk(ijkbeg:ijkend, 61) =  0.9999999999999999D-11
      rk(ijkbeg:ijkend, 62) =  DEXP(-0.2591722318817020D+02   &
         - ( -0.3310000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 63) =  DEXP(-0.3970958044115977D+02   &
        + ( 0.2000000000000000D+01 * LOG(temp(ijkbeg:ijkend)))   &
        +  0.9200000000000000D+02/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 64) =  0.1860000000000000D-11   &
        * DEXP(-( -0.1750000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.1320000000000000D-10
      rk(ijkbeg:ijkend, 65) =  DEXP(-0.2655601869289957D+02   &
         - ( -0.1900000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 66) =  DEXP(-0.2641908014195344D+02   &
         - ( -0.1900000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 67) =  0.1590000000000000D-13   &
        * DEXP(-( -0.5000000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.3800000000000000D-13
      rk(ijkbeg:ijkend, 68) =  DEXP(-0.2596142928067470D+02   &
         - (  0.2600000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 69) =  DEXP(-0.2870983077730048D+02   &
         - (  0.1900000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 70) =  DEXP(-0.2729454887930734D+02   &
         - (  0.1900000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 71) =   0.1620000000000000D-11   &
         * DEXP(-(  0.1900000000000000D+04 )/temp(ijkbeg:ijkend)) + & 
                 0.0000000000000000D+00   &
         * DEXP(-(  0.1500000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.1940000000000000D-13   &
         * DEXP(-(  0.1000000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 72) =  0.4920000000000000D-15
      rk(ijkbeg:ijkend, 73) =  0.4350000000000000D-17   &
         * DEXP( 0.2000000000000000D+01 * LOG(temp(ijkbeg:ijkend))   &
        -  0.2282000000000000D+04/temp(ijkbeg:ijkend)) + & 
                 0.1910000000000000D-13   &
         * DEXP(-(  0.4500000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.1080000000000000D-14   &
         * DEXP(-( -0.4500000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.0000000000000000D+00
      rk(ijkbeg:ijkend, 74) =  0.4000000000000000D-11   &
         * DEXP(-(  0.4460000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.0000000000000000D+00   &
         * DEXP(-( -0.4900000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 75) =  DEXP(-0.3551694253050293D+02   &
         - (  0.5000000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 76) =  0.8170000000000000D-14   &
         * DEXP(-(  0.2580000000000000D+04 )/temp(ijkbeg:ijkend)) + & 
                 0.4320000000000000D-15   &
         * DEXP(-(  0.1800000000000000D+04 )/temp(ijkbeg:ijkend)) + & 
                 0.2870000000000000D-16   &
         * DEXP(-(  0.8450000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.0000000000000000D+00   &
         * DEXP(-(  0.2283000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 77) =  0.7860000000000001D-14   &
         * DEXP(-(  0.1913000000000000D+04 )/temp(ijkbeg:ijkend)) + & 
                 0.0000000000000000D+00   &
         * DEXP(-(  0.7320000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 78) =  0.0000000000000000D+00   &
        * DEXP(-(  0.2112000000000000D+04 )/temp(ijkbeg:ijkend)) + & 
                 0.1380000000000000D-18
      rk(ijkbeg:ijkend, 79) =  DEXP(-0.3716986555487676D+02   &
         - (  0.1700000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 80) =  0.2000000000000000D-10
      rk(ijkbeg:ijkend, 81) =  0.9999999999999999D-11
      rk(ijkbeg:ijkend, 82) =  0.3600000000000000D-10
      rk(ijkbeg:ijkend, 83) =  DEXP(-0.3863712897853033D+02   &
         - ( -0.1044000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 83) = rk(ijkbeg:ijkend, 83) * SumM(ijkbeg:ijkend) * 0.2d0
      rk(ijkbeg:ijkend, 84) =  0.2800000000000000D-10
      Effko(ijkbeg:ijkend) =  0.9700000000000000D-28* (temp(ijkbeg:ijkend) / 3.d2)   &
                   **(- ( 0.5600000000000000D+01))
      Rapk(ijkbeg:ijkend) =  0.9300000000000000D-11* (temp(ijkbeg:ijkend) / 3.d2)   &
                    **(- ( 0.1500000000000000D+01))
      rk(ijkbeg:ijkend, 85) = (Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) /    &
                    ( 1.0d0 + Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend))) *   &
                    0.5860D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend)))**2))
      rk(ijkbeg:ijkend, 86) =  DEXP( 0.6462080260895155D+02   &
         - (  0.1395400000000000D+05 )/temp(ijkbeg:ijkend))
      Effko(ijkbeg:ijkend) =  0.9700000000000000D-28* (temp(ijkbeg:ijkend) / 3.d2)   &
                   **(- ( 0.5600000000000000D+01))
      Rapk(ijkbeg:ijkend) =  0.9300000000000000D-11* (temp(ijkbeg:ijkend) / 3.d2)   &
                    **(- ( 0.1500000000000000D+01))
      facteur(ijkbeg:ijkend) = (Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / ( 1.0d0 + Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) /    &
                 Rapk(ijkbeg:ijkend))) * 0.6d0 ** (1.0d0 / (1.0d0 +    &
                  (LOG10(Effko(ijkbeg:ijkend) * SumM(ijkbeg:ijkend) / Rapk(ijkbeg:ijkend)))**2))
      rk(ijkbeg:ijkend, 86) = facteur(ijkbeg:ijkend) * rk(ijkbeg:ijkend, 86)
      rk(ijkbeg:ijkend, 87) =  DEXP(-0.2619593659063922D+02   &
         - ( -0.1800000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 88) =  0.4360000000000000D-11
      rk(ijkbeg:ijkend, 89) =  0.6930000000000000D-11
      rk(ijkbeg:ijkend, 90) =  0.4000000000000000D-11
      rk(ijkbeg:ijkend, 91) =  0.4000000000000000D-11
      rk(ijkbeg:ijkend, 92) =  0.1220000000000000D-10
      rk(ijkbeg:ijkend, 93) =  0.4000000000000000D-11
      rk(ijkbeg:ijkend, 94) =  DEXP(-0.2859860514219026D+02   &
         - ( -0.8000000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 95) =  0.6160000000000000D-13   &
         * DEXP(-( -0.7000000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.1520000000000000D-12   &
         * DEXP(-( -0.1300000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 96) =  DEXP(-0.2934027936364486D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 98) =  DEXP(-0.2861185036894027D+02   &
         - ( -0.9800000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend, 99) =   0.5940000000000000D-12   &
         * DEXP(-( -0.5500000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.1990000000000000D-15   &
         * DEXP(-( -0.2640000000000000D+04 )/temp(ijkbeg:ijkend)) + & 
                 0.5560000000000000D-13   &
         * DEXP(-( -0.1300000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,100) =  DEXP(-0.2942678860655414D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,101) =  DEXP(-0.3002791688839384D+02   &
         - ( -0.4160000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,102) =  0.1030000000000000D-13   &
         * DEXP(-( -0.1580000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.6240000000000000D-13   &
         * DEXP(-( -0.4310000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.1530000000000000D-13   &
         * DEXP(-( -0.4670000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.4340000000000000D-14   &
         * DEXP(-( -0.6330000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,103) =  DEXP(-0.2948253058956238D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,104) =  DEXP(-0.2962612150917463D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,105) =  DEXP(-0.3096643075705270D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,106) =   0.1770000000000000D-10   &
         * DEXP(-(  0.4400000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.1480000000000000D-15   &
         * DEXP(-( -0.2510000000000000D+04 )/temp(ijkbeg:ijkend)) + & 
                 0.3100000000000000D-12   &
         * DEXP(-( -0.5080000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,107) =  DEXP(-0.2982027752361559D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,108) =  0.4440000000000000D-13   &
         * DEXP(-( -0.2110000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.2230000000000000D-12   &
         * DEXP(-( -0.4600000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.4100000000000000D-13   &
         * DEXP(-( -0.5220000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.1170000000000000D-13   &
         * DEXP(-( -0.6830000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,109) =  DEXP(-0.2846113415156165D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,110) =  DEXP(-0.2790545796163031D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,111) =  DEXP(-0.2864437356064584D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,112) =  0.7730000000000000D-12   &
         * DEXP(-( -0.5300000000000000D+03 )/temp(ijkbeg:ijkend)) + & 
                 0.1700000000000000D-12   &
         * DEXP(-( -0.5650000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,113) =  DEXP(-0.2835462750397320D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,114) =  DEXP(-0.3310607566097664D+02   &
         - ( -0.1000000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,115) =  DEXP(-0.3132793274173975D+02   &
         - ( -0.1000000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,116) =  0.1200000000000000D-11
      rk(ijkbeg:ijkend,117) =  0.1200000000000000D-11
      rk(ijkbeg:ijkend,118) =  0.1200000000000000D-11
      rk(ijkbeg:ijkend,119) =  0.1200000000000000D-11
      rk(ijkbeg:ijkend,120) =  0.1200000000000000D-11
      rk(ijkbeg:ijkend,121) =  0.3480000000000000D-11
      rk(ijkbeg:ijkend,122) =  0.1200000000000000D-11
      rk(ijkbeg:ijkend,123) =  DEXP(-0.2942678860655414D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,124) =  DEXP(-0.3274868498278332D+02   &
         - ( -0.1510000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,125) =  DEXP(-0.3171146277298166D+02   &
         - ( -0.1560000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,126) =  DEXP(-0.3717963534647257D+02   &
         - ( -0.2950000000000000D+04 )/temp(ijkbeg:ijkend))
      rk(ijkbeg:ijkend,127) =  0.4000000000000000D-11
      rk(ijkbeg:ijkend,128) =  0.1200000000000000D-11
    !   END DO
 
   END SUBROUTINE kinetic
 
  END MODULE mod_chem_spack_kinetic
 