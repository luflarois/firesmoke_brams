module modReadModisVeg

    use mem_grid, only: &
        nxtnest, nnxp, nnyp, npatch, nzp, nzg, ipm, jpm, platn, plonn, &
        iyear1, imonth1, idate1, ihour1, xtn, ytn, deltaxn,  deltayn, time, &
        deltax, deltay

    use node_mod, only: &
        mynum,         &
        mchnum,        &
        master_num

    use mem_mksfc, only: &
        dato, datp, cdato, nump, numpind, numpind1, numpind2, &
        ptable, glatp, glonp, datq_patch, npq

    use aer1_list, only: &
        nspecies

    use mem_aer1, only: aer1_g ! aer1_g(imode,ispc,1)%sc_p)
    !use mem_chem, only: chem1_g !chem1_g(ispc)%sc_p

    use node_mod, only: &
        nodei0, & ! intent(in)
        nodej0    ! intent(in)

    implicit none
    private

    logical, parameter :: vegDump =.false.
    character(len=256), parameter :: h5name = " "
    integer, parameter :: maxmiss=1000, ifm=1
    integer, parameter :: nveg=17, nveg_olson=97
    real, parameter :: fx = 1.0
    integer, parameter :: nspecies_aem=7, nveg_agreg=4
    integer, parameter :: nspecies_sfire = 6
    integer, parameter :: nfuelcats = 13
    character(LEN=4),dimension(nspecies_sfire),parameter :: sf_spc_name= &
    (/                     &
        'p25i', &
        'p25j', &
        'oc1 ', &
        'oc2 ', &
        'bc1 ', &
        'bc2 ', &    
    /) 	 
    
    integer, parameter :: p25i = 1
    integer, parameter :: p25j = 2
    integer, parameter :: oc1  = 3
    integer, parameter :: oc2  = 4
    integer, parameter :: bc1  = 5
    integer, parameter :: bc2  = 6

    character (len=29), parameter, dimension(nfuelcats+1) :: fuel_name = (/ &
            'Short grass (1 ft)           ', &
            'Timber (grass and understory)', &
            'Tall grass (2.5 ft)          ', &
            'Chaparral (6 ft)             ', &
            'Brush (2 ft)                 ', &
            'Dormant brush, hardwood slash', &
            'Southern rough               ', &
            'Closed timber litter         ', &
            'Hardwood litter              ', &
            'Timber (litter + understory) ', &
            'Light logging slash          ', &
            'Medium logging slash         ', &
            'Heavy logging slash          ', &
            'no fuel                      ' &
    /)
    

    integer, allocatable :: vegType(:,:,:,:)
    real,    dimension(nveg)   :: fc,bas_queiparms
    real :: burntArea

    type sfd
	    real    :: qarea 
        !# area
	    integer :: comb
        !# Comb type - From sfire
	    real    :: qtim
        !# burn time
	    real    :: frp
        !Fire radiate power
    end type sfd
    type(sfd), allocatable :: sfData(:,:,:)
    !# Data for fire point
    integer :: nSfirePoints
    !# Total of sfire points inside one grid point of BRAMS

    real,    dimension(nfuelcats,nspecies_sfire) :: emission_factor

    public :: readModisVeg, getEmission_per_pixel,nspecies_aem, AeM_spc_name, sfData

    !Emission factor from sfire
    data emission_factor/  &
    !--------------------------------------------------------------------------------------------------+
    !                              Emission Factor  by fuel cats                                       |
    !-------+------+------+------+------+------+------+------+------+------+------+------+------+------+
    !NFFL # |   1  |   2  |   3  |   4  |   5  |   6  |   7  |   8  |   9  |  10  |  11  |  12  |  13  |
    !-------+------+------+------+------+------+------+------+------+------+------+------+------+------+
               6.30, 13.00,  3.00, 21.00,  7.90,  5.00, 12.80,  7.00,  7.00,  7.00,  7.00,  7.00,  7.00, & !p25i
               6.30, 13.00,  3.00, 21.00,  7.90,  5.00, 12.80,  7.00,  7.00,  7.00,  7.00,  7.00,  7.00, & !p25j
               2.30,  1.30,  1.30, 11.00,  3.30,  3.30,  5.40,  4.60,  4.60,  4.60,  3.90,  3.90,  3.90, & !oc1
               2.30,  1.30,  1.30, 11.00,  3.30,  3.30,  5.40,  4.60,  4.60,  4.60,  3.90,  3.90,  3.90, & !oc2
               0.75,  0.30,  0.30,  0.40,  0.40,  0.40,  0.40,  0.45,  0.45,  0.45,  0.16,  0.16,  0.16, & !bc1
               0.75,  0.07,  0.07,  0.10,  0.10,  0.10,  0.10,  0.11,  0.11,  0.11,  0.04,  0.04,  0.04, & !bc2
    !-------+------+------+------+------+------+------+------+------+------+------+------+------+------+

contains

    subroutine readModisVeg(vegFile)


        character(len=*), intent(in) :: vegFile

        integer :: datq,lsp  &
        ,datsoil,soil_count,nlev_soil,nmiss      &
        ,jr,jp,ir,ip,iqv,ing,maxdq,jng,jng1,ngstor1    &
        ,ngstor2,npatpixs,nwat,ipat,idatq,isoil,jsoil,k      &
        ,no_veg ,iblksizo_veg ,isbego_veg ,iwbego_veg

        real :: deltallo_veg ,deltallo_soil ,deltallo_ndvi  &
        ,offlat_veg   ,offlat_soil   ,offlat_ndvi    &
        ,offlon_veg   ,offlon_soil   ,offlon_ndvi    &
        ,xp,yp,fracwat,plpp

        character(len=256) :: fnmiss(maxmiss)
        character(len=256) :: ivegtfn 

        ivegtfn = trim(vegFile)

        if (mchnum==master_num .and. time==0.0) then

            !Find size of patch arrays'
            call patch_array_size(npq, (xtn(2,1)-xtn(1,1)),            &
            1, ivegtfn, 0, 0, &
            0, 0)

            !Allocate arrays that need npq dimension',npq
            allocate (glatp(npq,npq,nnxp(ifm),nnyp(ifm)),  &
            glonp(npq,npq,nnxp(ifm),nnyp(ifm)),       &
            datq_patch(npq,npq,nnxp(ifm),nnyp(ifm)),  &
            datp(npq,npq,nnxp(ifm),nnyp(ifm)))
            
            !Allocate the array that will have the vegetation type
            allocate (vegType(npq,npq,nnxp(ifm),nnyp(ifm)))

            call patch_latlon(nnxp(ifm), nnyp(ifm),                  &
            xtn(1,ifm), ytn(1,ifm), deltaxn(ifm), deltayn(ifm), &
            platn(ifm), plonn(ifm))

            !print *,'LFR-DBG nnxp(1),nnyp(1)=',nnxp(1),nnyp(1),npq
            !print *, 'Reading header of MODIS vegetation types'
            call read_header(ivegtfn,iblksizo_veg,no_veg,isbego_veg  &
            ,iwbego_veg,offlat_veg,offlon_veg,deltallo_veg,'veg', h5name)

            !print *, 'Reading MODIS vegetation types in patchs'
            call fill_datp(nnxp(1),nnyp(1),no_veg,iblksizo_veg,isbego_veg,iwbego_veg  &
            ,platn,plonn,offlat_veg,offlon_veg,deltallo_veg,ivegtfn,iaction  &
            ,nmiss,fnmiss, h5name)

            vegType = datp

            if (vegDump) then
                do  ip = 1, npq
                    do jp = 1, npq
                        do ir = 1, nnxp(ifm)
                            do jr = 1, nnyp(ifm)
                                print *,ip,jp,ir,jr,vegType(ip,jp,ir,jr)
                            end do
                        end do
                    end do
                end do
            end if

            call setup_fireVeg()

        end if

    end subroutine readModisVeg


    subroutine getEmission_per_pixel(iPnt,jPnt,qsc)
        integer, intent(in) :: iPnt
        integer, intent(in) :: jPnt
        real, intent(out) :: qsc(nspecies_aem)

        integer :: iVeg
        integer :: ip,jp,ns, inpq, ia

        real :: qsc_aux(nspecies_aem,npq*npq)
        integer :: lIPnt, lJPnt

        inpq = 0
        qsc = 0
        lIPnt = nodei0(mynum,1)-iPnt
        lJPnt = nodej0(mynum,1)-jPnt

        !Loop over all points inside one BRAMS grid point
        do i=1,nSfirePoints
            call get_emission_per_One_point(i,qsc_aux(:,i))
        end do
        !Sum all emissions for each point inside a BRAMS grid point
        do i = 1,nSfirePoints
            do ns= 1,nspecies_aem
                qsc(ns) = qsc(ns) + qsc_aux(ns,i)
            end do
        end do

    end subroutine getEmission_per_pixel

    subroutine get_emission_per_One_point(iVeg,qsc)

        integer, intent(in) :: iVeg
        real,  intent(out), dimension(nspecies_aem) :: qsc
        integer :: i, iVeg_Agreg, iEsp
        real :: bas_by_fire, factArea, qarea, biomass_burned, qflam

        bas_by_fire = 0.0
        qsc = 0.0
        
        call agreg_veg(iVeg,iVeg_Agreg)
        if(iveg_agreg == 0) return

        bas_by_fire = bas_queiparms(iveg) !No Olson

        biomass_burned =  fc(iveg)*bas_by_fire*burntArea ! kg[dry biomass burned]	
        biomass_burned =  fx*biomass_burned !Adjust emissions
        qflam    =  flaming(iveg_agreg)

        !- avoid negative values for emission factor
        where(emission_factor <0.) emission_factor=0. 
        
        do iesp=1,nspecies_aem
            qsc(iesp) =  biomass_burned * emission_factor(iveg_agreg,iesp)*1.e-3 !kg[spc]
        end do

    end subroutine get_emission_per_One_point

    subroutine  agreg_veg(qveg,qveg_agreg)
        integer, intent(in)  :: qveg
        integer, intent(out) :: qveg_agreg

        integer :: catb(1:17)
        data catb/ &
               2, 1, 2, 1                 & !floresta tropical 2 and 4 / extra trop fores 1,3,5
             , 2, 3, 3, 3, 3              & !cerrado/woody savanna :6 a 9
             , 4, 4, 4, 4, 4, 0, 4, 0     / !pastagem/lavouras: 10 ...
        
        
        qveg_agreg = catb(qveg)

    end subroutine  agreg_veg


    subroutine setup_fireVeg()
        integer :: i
        real :: factArea, qarea

        do i=1,nveg
            fc(i)             = queiparms(1,i)        ! fracao de biomassa queimada
            bas_queiparms(i)  = queiparms(2,i)        ! kg/m2
        end do  

        allocate(sfData(nnxp(ifm),nnyp(ifm),nSfirePoints))
        sfcdata(:,:,:)%qarea  = 0.0
        sfcdata(:,:,:)%vtype = 0
        sfcdata(:,:,:)%qtim = 0.0
        sfcdata(:,:,:)%frp = 0.0

    end subroutine setup_fireVeg


end module modReadModisVeg