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
        'bc2 '  &    
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
    
    real,    dimension(nveg)   :: fc,bas_queiparms
    real :: burntArea
    integer :: nSfirePoints

    type sfd
	    real    :: FsizeTotal 
        !# area km2
	    integer :: comb
        !# Comb type - From sfire
	    real    :: qtim
        !# burn time
	    real    :: fre
        !Energia radiativa do fogo (Representa o total de biomassa
        !                             consumida pelo fogo). [MJ]
    end type sfd
    type(sfd), allocatable :: sfData(:,:,:)
    !# Data for fire point
 
    !# Total of sfire points inside one grid point of BRAMS

    real,    dimension(nfuelcats,nspecies_sfire) :: emission_factor

    public :: getEmission_per_pixel, sfData, setup_fire, sf_spc_name, nSfirePoints, nspecies_sfire

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
               0.75,  0.07,  0.07,  0.10,  0.10,  0.10,  0.10,  0.11,  0.11,  0.11,  0.04,  0.04,  0.04/   !bc2
    !-------+------+------+------+------+------+------+------+------+------+------+------+------+------+

contains

    subroutine getEmission_per_pixel(iPnt,jPnt,qsc,dn0,vol)
        integer, intent(in) :: iPnt
        integer, intent(in) :: jPnt
        real, intent(out) :: qsc(nspecies_sfire)
        real, intent(in) :: dn0
        real, intent(in) :: vol

        integer :: iVeg
        integer :: ip,jp,nf, inpq, ia

        real :: qsc_aux(nspecies_sfire,nSfirePoints)
        integer :: lIPnt, lJPnt, ns

        inpq = 0
        qsc = 0
        lIPnt = nodei0(mynum,1)-iPnt
        lJPnt = nodej0(mynum,1)-jPnt

        !Loop over all points inside one BRAMS grid point
        do ns=1,nSfirePoints
            !If there is no fire in this Sfire subgrid cycle 
            if(sfData(iPnt,jPnt,ns)%comb < 1 .or. sfData(iPnt,jPnt,ns)%comb > 13) cycle

            call get_emission_per_One_point(iPnt,jPnt,ns,qsc_aux(:,ns))
        end do
        !Sum all emissions for each point inside a BRAMS grid point
        do ns = 1,nSfirePoints
            !If there is no fire in this Sfire subgrid cycle 
            if(sfData(iPnt,jPnt,ns)%comb < 1 .or. sfData(iPnt,jPnt,ns)%comb > 13) cycle
            do nf= 1,nspecies_sfire
                !Divide the amount of src by air density, volume and time of emission
                qsc(nf) = qsc(nf) + qsc_aux(nf,ns)/dn0/vol/sfData(iPnt,jPnt,ns)%qtim
            end do
        end do     
        
    end subroutine getEmission_per_pixel

    subroutine get_emission_per_One_point(iPnt,jPnt,fPnt,qsc)

        integer, intent(in) :: iPnt,jPnt,fPnt
        real,  intent(out), dimension(nspecies_sfire) :: qsc
        integer :: nf, iVeg_Agreg, iEsp
        real :: bas_by_fire, factArea, qarea, biomass_burned, qflam

        ! FEER COEFICIENT for TPM
        biomass_burned =  sfData(iPnt,jPnt,fPnt)%fre! kg[dry biomass burned]
        
        qsc = 0.0

        !qflam(ifoc)    =  flaming(iveg_agreg(ifoc))

        !if(qfre > 0.) then
        !   qfires = 1.
        !else
        !   qfires = 0.
        !endif
        biomass_burned =  fx*biomass_burned

        do iesp=1,nspecies_sfire
            !!LFR - Retirei a divisão por emission_factor(iveg_agreg(ifoc),itpm) - verificar
            qsc(iesp) =  qsc(iesp) + biomass_burned * emission_factor(sfData(iPnt,jPnt,fPnt)%comb,iesp) !/emission_factor(iveg_agreg(ifoc),itpm)) !kg[spc]
        enddo

    end subroutine get_emission_per_One_point



    subroutine setup_fire()

        integer :: i,j,ns

        allocate(sfData(nnxp(ifm),nnyp(ifm),nSfirePoints))

        do i = 1, nnxp(ifm)
            do j = 1, nnyp(ifm)
                do ns = 1, nSfirePoints
                    sfData(i,j,ns)%FsizeTotal  = 0.0
                    sfData(i,j,ns)%comb = 0
                    sfData(i,j,ns)%qtim = 0.0
                    sfData(i,j,ns)%fre = 0.0
                end do
            end do 
        end do

    end subroutine setup_fire


end module modReadModisVeg