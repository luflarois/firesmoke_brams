!############################# Change Log ##################################
! Interface entre BRAMS e SFIRE
! 

SUBROUTINE sfclyr_sfire(mzp,mxp,myp,ia,iz,ja,jz)

     USE ModNamelistsfireFile, only : &
          namelistsfireFile, &
          CreateNamelistsfireFile, &
          DestroyNamelistsfireFile, &
          GetNamelistsfireFileName, &
          ReadNamelistsfireFile, &
          DumpNamelistsfireFile

     use node_mod, only:  &	
          mynum, &
          mchnum, &
          master_num, &
          nmachs, & !INTRUDOZIDO ISILDA
          nodemxp, & !INTRUDOZIDO ISILDA
          nodemyp, & !INTRUDOZIDO ISILDA
          nxbeg, & !INTRUDOZIDO ISILDA
          nxend, & !INTRUDOZIDO ISILDA
          nybeg, & !INTRUDOZIDO ISILDA
          nyend !INTRUDOZIDO ISILDA     
     USE mem_grid ,  ONLY : ngrid, deltaxn, deltayn, dtlongn,time,istp, &
                            nnzp,nnxp,nnyp,grid_g,npatch,zt,nzpmax           ! dtlong, ESTES SO INTROD POR ISILDA

     USE mem_sfire, only :   &
          sfire_g,     &
          alloc_sfire_brams, &
          nullify_sfire_brams, &
          dealloc_sfire_brams, &
          zero_sfire_brams, &
          config_flags

     use module_domain_type
     use module_fr_sfire_driver_brams
     use module_model_constants    	
     USE, INTRINSIC :: IEEE_ARITHMETIC
     USE ReadBcst, ONLY: &
            gatherData,Broadcast
     use ParLib, only:  parf_barrier


     USE mem_cuparm, only: cuparm_g 
     USE mem_stilt, only: stilt_g
     USE mem_micro, only: micro_g
     USE mem_basic, only: basic_g
     USE mem_leaf, only: leaf_g
     USE mem_jules, only: jules_g
     USE mem_varinit, only: varinit_g  
     USE mem_turb,  only : turb_g
     USE rconstants , only: cpi,p00,cpor
     USE mem_basic   , only: basic_g
     USE micphys, only: mcphys_type     
     use modreadModisVeg, only: readModisVeg, getEmission_per_pixel

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: mzp,mxp,myp,ia,iz,ja,jz
     ! REAL, INTENT(IN) :: time
     REAL :: fdx, fdy

     !type(grid_config_rec_type), pointer :: config_flags => null()  

     !Local Variables
     INTEGER            :: ng,a_step,i ,j, k, ierr,np 	
     real ::hf(nnxp(1),nnyp(1))
     real ::ch(nnxp(1),nnyp(1))
     real ::qf(nnxp(1),nnyp(1))
     real ::cq(nnxp(1),nnyp(1)) 
     real ::temp_sflux_t(nnxp(1),nnyp(1))
     real ::temp_sflux_r(nnxp(1),nnyp(1))
     ! real ::temp_aconpr(nnxp(1),nnyp(1))
     real ::temp_accpr(nnxp(1),nnyp(1))
     real ::temp_accpp(nnxp(1),nnyp(1))
     real ::temp_accps(nnxp(1),nnyp(1))
     real ::temp_accpa(nnxp(1),nnyp(1))
     real ::temp_accpg(nnxp(1),nnyp(1))
     real ::temp_accph(nnxp(1),nnyp(1))
     real ::temp_t2mj(nnxp(1),nnyp(1))
     real ::temp_topt(nnxp(1),nnyp(1))
     real ::temp_rv2mj(nnzp(1),nnxp(1),nnyp(1))
     real ::temp_glat(nnxp(1),nnyp(1))
     real ::temp_glon(nnxp(1),nnyp(1))
     real ::temp_veg_rough(nnxp(1),nnyp(1),npatch)
     real ::temp_veg_class(nnxp(1),nnyp(1),npatch)
     ! real ::zt(nzpmax)
     real ::temp_rtgt(nnxp(1),nnyp(1))
     real ::temp_h(nnzp(1),nnxp(1),nnyp(1))
     !  real ::temp_dnp(nnzp(1),nnxp(1),nnyp(1))
     real ::temp_up(nnzp(1),nnxp(1),nnyp(1))
     real ::temp_vp(nnzp(1),nnxp(1),nnyp(1))
     real ::temp_dn0(nnzp(1),nnxp(1),nnyp(1))
     CHARACTER(len=16)  :: varn
     real :: reduz
     real ::temp_pi0(nnzp(1),nnxp(1),nnyp(1))
     real ::temp_pp(nnzp(1),nnxp(1),nnyp(1))
     real ::temp_press(nnxp(1),nnyp(1))
     real ::temp_theta(nnzp(1),nnxp(1),nnyp(1))
     real ::picpi
     real ::a,b,step_isil
     integer, parameter :: ifm = 1 

     varn = 'SFLUX_T'
     call gatherData(2, varn,1, nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          turb_g(ngrid)%sflux_t, temp_sflux_t)
     varn = 'SFLUX_R'
     call gatherData(2, varn,1, nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          turb_g(ngrid)%sflux_r, temp_sflux_r)
     varn = 'ACCPR'
     call gatherData(2, varn,1, nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          micro_g(ngrid)%accpr, temp_accpr)
     varn = 'ACCPP'
     call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          micro_g(ngrid)%accpp, temp_accpp)
     varn = 'ACCPS'
     call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          micro_g(ngrid)%accps, temp_accps)
     varn = 'ACCPA'
     call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          micro_g(ngrid)%accpa, temp_accpa)
     varn = 'ACCPG'
     call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          micro_g(ngrid)%accpg, temp_accpg)
     varn = 'ACCPH'
     call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          micro_g(ngrid)%accph, temp_accph)
     varn = 'TOPT'
     call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          grid_g(ngrid)%topt, temp_topt)
     varn = 'RV'  !***ATENCAO acoplamos o leaf
     call gatherData(3, varn, 1, nnzp(1), nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          basic_g(ngrid)%rv, temp_rv2mj)	    
     varn = 'PP'
     call gatherData(3, varn, 1, nnzp(1),nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          basic_g(ngrid)%pp, temp_pp)
     varn = 'GLAT'
     call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          grid_g(ngrid)%glat, temp_glat)
     varn = 'GLON'
     call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          grid_g(ngrid)%glon, temp_glon)
     do np=1,npatch
          varn = 'VEG_ROUGH'
          call gatherData(2, varn, 1,nnxp(1), nnyp(1), &
            nmachs, mchnum, mynum, master_num,                    &
            leaf_g(ngrid)%veg_rough(1:nnxp(1),1:nnyp(1),np), &
	    temp_veg_rough(1:nnxp(1),1:nnyp(1),np))
     enddo  
     varn = 'LEAF_CLASS'
     do np=1,npatch
          call gatherData(2, varn, 1,nnxp(1), nnyp(1), &
            nmachs, mchnum, mynum, master_num, &
            leaf_g(ngrid)%leaf_class(1:nnxp(1),1:nnyp(1),np), &
	       temp_veg_class(1:nnxp(1),1:nnyp(1),np))
     enddo   
     varn = 'RTGT'
     call gatherData(2, varn, 1, nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          grid_g(ngrid)%rtgt, temp_rtgt)	    
     varn = 'UP'
     call gatherData(3, varn, 1, nnzp(1),nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          basic_g(ngrid)%up, temp_up)
     varn = 'VP'
     call gatherData(3, varn, 1, nnzp(1),nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          basic_g(ngrid)%vp, temp_vp)
     varn = 'DN0'
     call gatherData(3, varn, 1, nnzp(1),nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          basic_g(ngrid)%dn0, temp_dn0)
     varn = 'THETA'
     call gatherData(3, varn, 1, nnzp(1),nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          basic_g(ngrid)%theta, temp_theta)
     varn = 'PI0'
     call gatherData(3, varn, 1, nnzp(1),nnxp(1), nnyp(1), &
          nmachs, mchnum, mynum, master_num,                    &
          basic_g(ngrid)%pi0, temp_pi0)


  if (mchnum==master_num) then

     if (time==0.0) then
          call readModisVeg('/mnt/beegfs/luiz.rodrigues/bin/datafix_model/veg_modis/MODIS') !LFR DBG
     endif
     if (time==30.0) then
      call getEmission_per_pixel(45,90)
     end if

  !***********************************************************************************
  !***********************************************************************************
  !***********************************************************************************
  goto 1000
   !***********************************************************************************
  !***********************************************************************************
  !*********************************************************************************** 


  !print*,'GLAT',temp_glat(1,1),temp_glat(1,200), temp_glat(200,1),temp_glat(200,200)
  !print*,'GLON',temp_glon(1,1),temp_glon(1,200), temp_glon(200,1),temp_glon(200,200)

!  print*,'LIX0 1 *********'
 !      print*,temp_veg_class(:,:,1)
  !     print*,'LIXO 2'
   !    print*,temp_veg_class(:,:,2)
    !    print*,'LIXO 3'
     !  print*,temp_veg_class(:,:,3)
     !    print*,'LIXO 4'
     !  print*,temp_veg_class(:,:,4)
  do k=1,nnzp(1)
    do j=1,nnyp(1)
      do i=1,nnxp(1)
         temp_h(k,i,j)=zt(k) * temp_rtgt(i,j)
      enddo
    enddo
  enddo
  


 k=2
 !MFF - Troquei nnxp e nnyp de lugar
 do j=1,nnyp(1)
   do i=1,nnxp(1)
    picpi = (temp_pi0(k,i,j) + temp_pp(k,i,j)) * cpi
    temp_press(i,j) = p00 * picpi ** cpor
    temp_t2mj(i,j) = temp_theta(k,i,j) * picpi  !ACOPLAMENTO DO LEAF
   enddo
 enddo   


! temp_topt(:,:)=60.0 
 	 
  ! print*, " Iniciando driver do sfire " 
   call flush(6)
  ! print*,'DEBUG :: -----------------'  
  ! print*, 'JULES INICIO'
  !print*,'============='
  !print*,'DEBUG ::::temp_h::::'
  !print*,temp_h
  !print*,'============='
  !print*,'DEBUG ::::temp_veg_rough::::'
  !print*,temp_veg_rough
  !print*,'FIM DEBUG :: -----------------' 
   !--- Inicializando a grade do JULES ---
 
   print *, "DEBUG - ISTP ", istp
   IF((istp-1) == 0) THEN
  
   print *, "DEBUG - FIST ISTP ", istp
      call GlobalLatLonGrid(temp_glon,temp_glat)
     
      print*, " Criando namelist do SFIRE" 
      call flush(6)
      call CreateNamelistsfireFile(config_flags)
      
      print*, " Pegando namelist do SFIRE" 
      call flush(6)
      call GetNamelistsfireFileName(config_flags)
      
      
      print*, " Lendo namelist do SFIRE" 
      call flush(6)
      call ReadNamelistsfireFile(config_flags)
 
      WRITE(*,"(50('-'),/,a)")'Reading model control file (./sfire.in)...'
      
      
      
      ! get fire mesh dimensions
      CALL get_ijk_from_subgrid (  config_flags ) 
      
      
      sfire_g%sr_x = config_flags%sr_x
      sfire_g%sr_y = config_flags%sr_y
      
      sfire_g%dx=deltaxn(ngrid)
      sfire_g%dy=deltayn(ngrid)
      
      sfire_g%itimestep = istp-1
      sfire_g%num_tiles = 1
      sfire_g%num_tiles_spec = 1
   

   
     
      call nullify_sfire_brams(sfire_g)
      call alloc_sfire_brams(sfire_g,config_flags%ims,config_flags%ime, config_flags%kms,config_flags%kme, config_flags%jms,config_flags%jme &
		,config_flags%ifms,config_flags%ifme, config_flags%jfms,config_flags%jfme, config_flags%nfmc)
          ! Putting zero on all values
      call zero_sfire_brams(sfire_g)
      
      
      
      IF ( ASSOCIATED(sfire_g%i_start) ) THEN ; DEALLOCATE( sfire_g%i_start ) ; NULLIFY( sfire_g%i_start ) ; ENDIF
      IF ( ASSOCIATED(sfire_g%i_end) ) THEN ; DEALLOCATE( sfire_g%i_end ) ; NULLIFY( sfire_g%i_end ) ; ENDIF
      IF ( ASSOCIATED(sfire_g%j_start) ) THEN ; DEALLOCATE( sfire_g%j_start ) ; NULLIFY( sfire_g%j_start ) ; ENDIF
      IF ( ASSOCIATED(sfire_g%j_end) ) THEN ; DEALLOCATE( sfire_g%j_end ) ; NULLIFY( sfire_g%j_end ) ; ENDIF
      ALLOCATE(sfire_g%i_start(sfire_g%num_tiles))
      ALLOCATE(sfire_g%i_end(sfire_g%num_tiles))
      ALLOCATE(sfire_g%j_start(sfire_g%num_tiles))
      ALLOCATE(sfire_g%j_end(sfire_g%num_tiles))  
      
      sfire_g%i_start(1) = config_flags%ips
      sfire_g%i_end(1)   = config_flags%ipe
      sfire_g%j_start(1) = config_flags%jps
      sfire_g%j_end(1)   = config_flags%jpe
      
      print*, " Lendo dados de combustiveis do SFIRE -sim 0" 
      call flush(6)

      call combinit_user(temp_veg_class, temp_glat, temp_glon,sfire_g%nfuel_cat)
      
     ! print*,"sfire_g%nfuel_cat - sfcsdrive - no init"
     ! print*,sfire_g%nfuel_cat
      

      print*, " Convertendo indices das variaveis para o SFIRE - simu 0" 
      call flush(6)
      
      call swap_brams_sfire( sfire_g,& !temp_aconpr,&
           temp_accpr,temp_accpp,temp_accps,temp_accpa,&
           temp_accpg,temp_accph,temp_t2mj,temp_topt,&
           temp_rv2mj,temp_pp,temp_glat,temp_glon,&
           temp_veg_rough,temp_h,temp_dn0,temp_up,&
           temp_vp,temp_press )
     
    !gradiente da topografia (Atencao ao loop deve ser jfts no acplamento final)
     a=sfire_g%dx/sfire_g%sr_x
     b=sfire_g%dy/sfire_g%sr_y
     do i=2,nnxp(1)-1
      do j=2,nnyp(1)-1
       sfire_g%dzdxf(i,j)=(sfire_g%zsf(i+1,j)-sfire_g%zsf(i-1,j))/(2.*a)
       sfire_g%dzdyf(i,j)=(sfire_g%zsf(i,j+1)-sfire_g%zsf(i,j-1))/(2.*b)
      enddo
    enddo

   ! print*,'LIXO DRIVE1'
   ! print*,'DZDXF'
   ! print*, sfire_g%dzdxf
   ! print*,'DZDYF'
   ! print*,sfire_g%dzdyf


      print*, " Iniciando o SFIRE" 
      call flush(6)
      
      step_isil=(dtlongn(ngrid)*((config_flags%dx/1000.)*6.))/dtlongn(ngrid)
      
      print*, 'ISILDA1 INIT :: time, dtlongn', time, dtlongn
      call flush(6)
   
      call sfire_driver_em_init ( sfire_g , config_flags, time, step_isil )
      
      
      
      
      
      where( ieee_is_nan(sfire_g%rho)) sfire_g%rho = 0.00001
      where( sfire_g%rho <= 0.00000000000000000000001) sfire_g%rho = 0.00001
      
      where( ieee_is_nan(sfire_g%dz8w)) sfire_g%dz8w = 1.
      where( sfire_g%dz8w < 1.) sfire_g%dz8w = 1.
	    
   ENDIF !(time==0)
 
   sfire_g%itimestep = istp
   print*,"sfire_g%itimestep =",sfire_g%itimestep
   print*, " Convertendo indices das variaveis para o SFIRE - simu run" 
   call flush(6)
 
  ! get fire mesh dimensions
   CALL get_ijk_from_subgrid (  config_flags ) 
   
    print*, " sai do get_ijk"
    call flush(6)   
      
      sfire_g%sr_x = config_flags%sr_x
      sfire_g%sr_y = config_flags%sr_y
      
      sfire_g%dx=deltaxn(ngrid)
      sfire_g%dy=deltayn(ngrid)

      sfire_g%i_start(1) = config_flags%ips
      sfire_g%i_end(1)   = config_flags%ipe
      sfire_g%j_start(1) = config_flags%jps
      sfire_g%j_end(1)   = config_flags%jpe
  
   print*, " vou para o swap_brams_sfire"
   call flush(6)
   call swap_brams_sfire( sfire_g,&   !temp_aconpr,&
           temp_accpr,temp_accpp,temp_accps,temp_accpa,&
           temp_accpg,temp_accph,temp_t2mj,temp_topt,&
           temp_rv2mj,temp_pp,temp_glat,temp_glon,&
           temp_veg_rough,temp_h,temp_dn0,temp_up,&
           temp_vp,temp_press )
     
      print*, " sai do swap_brams_sfire"
   call flush(6)
      	   
     !gradiente da topografia (Atencao ao loop deve ser jfts no acplamento final)
     a=sfire_g%dx/sfire_g%sr_x
     b=sfire_g%dy/sfire_g%sr_y
     do i=2,nnxp(1)-1
      do j=2,nnyp(1)-1
       sfire_g%dzdxf(i,j)=(sfire_g%zsf(i+1,j)-sfire_g%zsf(i-1,j))/(2.*a)
       sfire_g%dzdyf(i,j)=(sfire_g%zsf(i,j+1)-sfire_g%zsf(i,j-1))/(2.*b)
      enddo
    enddo
	print*, " passei o calculo grad topografia"
   call flush(6)

   ! print*,'LIXO DRIVE2'
   ! print*,'DZDXF'
   ! print*, sfire_g%dzdxf
   ! print*,'DZDYF'
   ! print*,sfire_g%dzdyf

    sfire_g%rainc=0.
       
      where( ieee_is_nan(sfire_g%rho)) sfire_g%rho = 0.00001
      where( sfire_g%rho <= 0.00000000000000000000001) sfire_g%rho = 0.00001
      
      where( ieee_is_nan(sfire_g%dz8w)) sfire_g%dz8w = 1.
      where( sfire_g%dz8w < 1.) sfire_g%dz8w = 1.	   

   print*, " Integrando o SFIRE " 
   call flush(6)
   

   ! Comentando apenas para testar a inicializacao do modulo e compilacao
   step_isil=(dtlongn(ngrid)*((config_flags%dx/1000.)*6.))/dtlongn(ngrid)
   
   call sfire_driver_em_step (  sfire_g , config_flags, time, step_isil )
   
   !passagem das unidades W/m^2 para K m/s
  
      
    hf(:,:)=0.0
    ch(:,:)=0.0
    do j=config_flags%jps,config_flags%jpe
            do i=config_flags%ips,config_flags%ipe
           ! reduz = dn0(2,i,j) * 1004.
             reduz = sfire_g%rho(i,2,j) * 1004
            hf(i,j) = sfire_g%grnhfx(i,j) / reduz
            ch(i,j) = sfire_g%canhfx(i,j) / reduz
	    !if(hf(i,j) > 0.)then
	    !print*,"ISILDA4 - h"
	    !print*,hf(i,j)
	    !endif 
            enddo
    enddo
  
    
   ! print*,'isilda1-',time
   ! print*,'heat1'
   ! print*,sfire_g%grnhfx
   ! print*,'heat2'
   ! print*,hf
   
   !passagem das unidades de W/m^2 para Kg/kg m/s 
   
    qf(:,:)=0.0
    cq(:,:)=0.0
   do j=config_flags%jps,config_flags%jpe
            do i=config_flags%ips,config_flags%ipe
            !reduz = dn0(2,i,j) * 2.5e6
            reduz = sfire_g%rho(i,2,j) * 2.5e6
            qf(i,j) = sfire_g%grnqfx(i,j) / reduz
            cq(i,j) = sfire_g%canqfx(i,j) / reduz
	   ! if(qf(i,j) > 0.)then
	   ! print*,"ISILDA4 - le"
	   ! print*,qf(i,j)
	   ! endif 
            enddo
    enddo

    
  ! soma ao calor do  Jules


    do j=1,config_flags%jps - 1
            do i=1,nnxp(1)
   
    temp_sflux_t(i,j) = temp_sflux_t(i,j)
    temp_sflux_r(i,j) = temp_sflux_r(i,j)
            enddo
    enddo
    do j=config_flags%jpe + 1,nnyp(1)
            do i=1,nnxp(1)
   
    temp_sflux_t(i,j) = temp_sflux_t(i,j)
    temp_sflux_r(i,j) = temp_sflux_r(i,j)
            enddo
    enddo
    do j=config_flags%jps,config_flags%jpe 
            do i=1,config_flags%ips-1
   
    temp_sflux_t(i,j) = temp_sflux_t(i,j)
    temp_sflux_r(i,j) = temp_sflux_r(i,j)
            enddo
    enddo
    do j=config_flags%jps,config_flags%jpe 
            do i=config_flags%ipe+1,nnxp(1)
   
    temp_sflux_t(i,j) = temp_sflux_t(i,j)
    temp_sflux_r(i,j) = temp_sflux_r(i,j)
            enddo
    enddo
    do j=config_flags%jps,config_flags%jpe
            do i=config_flags%ips,config_flags%ipe
               temp_sflux_t(i,j) = temp_sflux_t(i,j) + hf(i,j) + ch(i,j)
             !  print*,"temp_sflux_t",temp_sflux_t(i,j)
             !  call flush(6) 
               temp_sflux_r(i,j) = temp_sflux_r(i,j) + qf(i,j) + cq(i,j)
             !  print*,"temp_sflux_r",temp_sflux_r(i,j)
             !  call flush(6)
             
            enddo
     enddo
  ! print*,'DEBUG :: -----------------'  
  ! print*,'jules'
  ! print*,temp_sflux_t
  ! print*,'FIM DEBUG :: ------------'
  !paraleliza
   !print*,"ISILDA"
   !print*,"TEMP_SFLUX_T"
   !print*,temp_sflux_t
   !print*,"ISILDA"
   !print*,"sfire_g%grnhfx"
   !print*,sfire_g%grnhfx
   !print*,"ISILDA"
   !print*,"TEMP_SFLUX_R"
   !print*,temp_sflux_r
   !print*,"ISILDA"
   !print*,"sfire_g%grnqfx"
   !print*,sfire_g%grnqfx 
  !***********************************************************************************
  !***********************************************************************************
  !***********************************************************************************
1000 continue
  !***********************************************************************************
  !***********************************************************************************
  !***********************************************************************************   
  endif
 ! print*,"estou no sfire vou entrar no PARALELISMO - ISTP - TIME =",istp,time 
!   call parf_barrier(0)
  ! print*,"dentro drive vou para Broadcast- temp_sflux_t"
!!!!!!   call Broadcast(temp_sflux_t, master_num,   'temp_sflux_t' )
  ! print*,"dentro drive sai do Broadcast- temp_sflux_t",mynum
  ! print*,"dentro drive vou para Broadcast- temp_sflux_R"
!!!!!   call Broadcast(temp_sflux_r, master_num,   'temp_sflux_r' )
  ! print*,"dentro drive sai do Broadcast- temp_sflux_R",mynum
  ! print*,"dentro drive vou para o swap"
  !call swap_sfire_brams( temp_sflux_t,temp_sflux_r, &
   !        turb_g(ngrid)%sflux_t,turb_g(ngrid)%sflux_r)
   !print*,"dentro drive sai do swap" 
  
      ! varn = 'SFLUX_T'
      ! call gatherData(2, varn,1, nnxp(1), nnyp(1), &
      !      nmachs, mchnum, mynum, master_num,                    &
      !       temp_sflux_t, turb_g(ngrid)%sflux_t)

        !    print*,"passei SFLUX_T"

      ! varn = 'SFLUX_R'
      ! call gatherData(2, varn,1, nnxp(1), nnyp(1), &
      !      nmachs, mchnum, mynum, master_num,                    &
      !      temp_sflux_r,turb_g(ngrid)%sflux_r)


     !  call swap_sfire_brams( temp_sflux_t,temp_sflux_r, &
     !      turb_g(ngrid)%sflux_t,turb_g(ngrid)%sflux_r)
      !-- scattering local data
       !print*, "turb_g(ngrid)%sflux_t", turb_g(ngrid)%sflux_t
       !call flush(6)
       !print*,"turb_g(ngrid)%sflux_r", turb_g(ngrid)%sflux_r
       !call flush(6)
    !  call mk_4_buff(turb_g(ngrid)%sflux_t, temp_sflux_t, &
    !        mzg, nnxp(ifm), nnyp(ifm), npat, mzg, n2, n3, npat, ia, iz, ja, jz)

     ! call mk_4_buff(turb_g(ngrid)%sflux_r, temp_sflux_r, &
     !       mzg, nnxp(ifm), nnyp(ifm), npat, mzg, n2, n3, npat, ia, iz, ja, jz)

   !   call mk_2_buff(turb_g(ngrid)%sflux_t, temp_sflux_t, &
    !                          nnxp(ifm), nnyp(ifm), &
    !                          nodemxp(mynum,ifm), nodemyp(mynum,ifm), &
    !                          nxbeg(mynum,ifm),nxend(mynum,ifm),nybeg(mynum,ifm),nyend(mynum,ifm))

    !  call mk_2_buff(turb_g(ngrid)%sflux_r, temp_sflux_r, &
     !                         nnxp(ifm), nnyp(ifm), &
     !                         nodemxp(mynum,ifm), nodemyp(mynum,ifm), &
      !                        nxbeg(mynum,ifm),nxend(mynum,ifm),nybeg(mynum,ifm),nyend(mynum,ifm))




 
  ! call StoreOwnChunk_2D(ngrid,temp_sflux_t, turb_g(ngrid)%sflux_t, nnxp(1), nnyp(1), 'SFLUX_T')	
  ! call StoreOwnChunk_2D(ngrid,temp_sflux_r, turb_g(ngrid)%sflux_r, nnxp(1), nnyp(1), 'SFLUX_R')
  ! print*,"DEPOIS DA introducao do PARALELISMO - vou sair do sfire"
  ! print*,"turb_g(ngrid)%sflux_t"
 !  print*,turb_g(ngrid)%sflux_t
  ! print*,"turb_g(ngrid)%sflux_r"
 !  print*,turb_g(ngrid)%sflux_r
   
 !  print*,"VOU SAIR DO SFIRE"  
   RETURN
   
END SUBROUTINE sfclyr_sfire

SUBROUTINE swap_brams_sfire( sfire,&      !temp_aconpr,&
           temp_accpr,temp_accpp,temp_accps,temp_accpa,&
           temp_accpg,temp_accph,temp_t2mj,temp_topt,&
           temp_rv2mj,temp_pp,temp_glat,temp_glon,&


 
           temp_veg_rough,temp_h,temp_dn0,temp_up,&
           temp_vp,temp_press)


  USE module_domain_type
  USE mem_grid ,  ONLY : nnzp,nnxp,nnyp, &
                      polelat,polelon,npatch

  implicit none
 
  !real, INTENT(IN) ::temp_aconpr(nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_accpr(nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_accpp(nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_accps(nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_accpa(nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_accpg(nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_accph(nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_t2mj(nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_topt(nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_rv2mj(nnzp(1),nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_pp(nnzp(1),nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_glat(nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_glon(nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_veg_rough(nnxp(1),nnyp(1),npatch)
  real, INTENT(IN) ::temp_h(nnzp(1),nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_dn0(nnzp(1),nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_up(nnzp(1),nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_vp(nnzp(1),nnxp(1),nnyp(1))
  real, INTENT(IN) ::temp_press(nnxp(1),nnyp(1))
  real :: x, y  
  


  TYPE(domain) , TARGET :: sfire  
  INTEGER :: k
  REAL, DIMENSION(nnzp(1),nnxp(1),nnyp(1)) :: temp
 

 !print*,"entrei no swap_brams_sfire"
 !print*,"rainc"
 ! sfire%rainc = max(temp_aconpr(:,:),0.0)
  sfire%rainc(:,:) = 0.0
 ! print*,"rainnc"
  sfire%rainnc(:,:) = sfire%rainc(:,:) + temp_accpr(:,:) &
  	+ temp_accpp(:,:) + temp_accps(:,:) &
	+ temp_accpa(:,:) + temp_accpg(:,:) &
	+ temp_accph(:,:)

 ! print*,"***"

  ! TO DO: Rever metodo para preencher variavel
  
 ! print*,"PASSOU RAINC E RAINNC"
 ! print*,"vou para o t2"
  sfire%t2(:,:) =  temp_t2mj(:,:)
  
  !sfire%t2(:,:) =  jules_g(ngrid)%t2mj(:,:)
  !print*,"vou para o ht"
  sfire%ht(:,:) = temp_topt(:,:)
  !print*,"vou para o zsf"
  sfire%zsf(:,:) = temp_topt(:,:)
  !print*,"vou para o q2"
  sfire%q2(:,:) = temp_rv2mj(2,:,:)
  !print*,"vou para o mut"
  sfire%mut(:,:) = temp_pp(1,:,:)
  !print*,"vou para o xlat"
  sfire%xlat(:,:) = temp_glat(:,:)
  !print*,"vou para o xlong"
  sfire%xlong(:,:) = temp_glon(:,:)  
  !print*,"vou para o fxlat"
  sfire%fxlat(:,:) = temp_glat(:,:)
  !print*,"vou para o fxlong"
  sfire%fxlong(:,:) = temp_glon(:,:)
  !print*,"vou para o z0"
  sfire%z0(:,:) = temp_veg_rough(:,:,2)
  sfire%psfc(:,:) = temp_press(:,:)
  !print*,"vou para a temp"
   do k=1, nnzp(1)-1
   temp(k,:,:) = temp_h(k+1,:,:) - temp_h(k,:,:)
  enddo
 ! print*,"TEMP - DZ8W"
  !print*,temp
  ! Isso mudar quando tiver topografia em alta resoluo do modelo de fogo (sr_x, sr_y)
  !sfire%dzdxf(:,:) = grid_g(ngrid)%glon(:,:)
  !sfire%dzdyf(:,:) = grid_g(ngrid)%glat(:,:)
  
  
  ! TO DO: Variaveis que precisam trocar indices
 ! print*,"vou para o call rho"
  call swap_kij_to_ikj(sfire%rho, temp_dn0)
  !print*,"vou para o call dz8w"
  call swap_kij_to_ikj(sfire%dz8w, temp)
  !print*,"vou para o call z_at_w"
  call swap_kij_to_ikj(sfire%z_at_w, temp_h)
  
  sfire%ph_2(:,:,:) = sfire%z_at_w(:,:,:)
  sfire%phb(:,:,:) = sfire%z_at_w(:,:,:)
  !print*,"vou para o call u_2"
  call swap_kij_to_ikj(sfire%u_2, temp_up)
  !print*,"vou para o call v_2"
  call swap_kij_to_ikj(sfire%v_2, temp_vp)
 
   !print*,"fiz tudo vou sair da sawp_fire"
   return
END SUBROUTINE swap_brams_sfire

 SUBROUTINE swap_sfire_brams( temp_sflux_t, temp_sflux_r, &
            said_sflux_t, said_sflux_r)


  use mem_grid, only: nnxp, nnyp    

  use node_mod, only:  &	
       mynum, &
       master_num, &
       nodemxp, &
       nodemyp, &
       ia, iz, &
       ja, jz
  USE ParLib, ONLY: &
         parf_bcast


  implicit none
  real, INTENT(IN) ::temp_sflux_t(nnxp(1), nnyp(1))
  real, INTENT(IN) ::temp_sflux_r(nnxp(1), nnyp(1))
 ! real, INtENT(OUT) ::said_sflux_t(nnxp(1), nnyp(1))
  real, INtENT(OUT) ::said_sflux_t(nodemxp(mynum,1), nodemyp(mynum,1))
 ! real, INTENT(OUT) ::said_sflux_r(nnxp(1), nnyp(1))
  real, INTENT(OUT) ::said_sflux_r(nodemxp(mynum,1), nodemyp(mynum,1))
  integer:: n1,n2

   n1=nnxp(1)
   n2=nnyp(1)
 
   !print*,"n1,n2,ia, iz, ja, jz =", n1,n2,ia, iz, ja, jz, nodemxp(mynum,1), nodemyp(mynum,1)
   !CALL parf_bcast(temp_sflux_t, n1, &
   !         n2, master_num)
    !print*,"dentro do swap - antes do mkbuff sflux_t"
    call flush(6)
       ! Distributing local information about sensible heat
   CALL mk_2_buff(temp_sflux_t(:,:), said_sflux_t(:,:), &
            n1, n2, nodemxp(mynum,1), nodemyp(mynum,1), &
            ia, iz, ja, jz)
   !print*,"ia, iz, ja, jz =",ia, iz, ja, jz	    
   !print*,"dentro do swap - depois do mkbuff sflux_t"
   call flush(6)
  ! CALL parf_bcast(temp_sflux_r, n1, &
   !         n2, master_num)
    !print*,"dentro do swap - antes do mkbuff sflux_r"
    call flush(6)
       ! Distributing local information about latent heat
   CALL mk_2_buff(temp_sflux_r(:,:), said_sflux_r(:,:), &
             n1, n2, nodemxp(mynum,1), nodemyp(mynum,1),&
             ia, iz, ja, jz)
   !print*,"ia, iz, ja, jz =",ia, iz, ja, jz	    
   !print*,"dentro do swap - depois do mkbuff sflux_r"
   call flush(6)
  return

END SUBROUTINE swap_sfire_brams

SUBROUTINE swap_kij_to_ikj(a, b)
 ! USE node_mod, only: mxp, myp, mzp
  USE mem_grid ,  ONLY :nnzp,nnxp,nnyp,ngrid
  implicit none
  
  integer :: i, j, k
  real :: b(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)), a(nnxp(ngrid),nnzp(ngrid),nnyp(ngrid))
  
  
  do j=1, nnyp(ngrid)
    do k=1, nnzp(ngrid)
      do i=1, nnxp(ngrid)
        a(i,k,j) = b(k,i,j)
      enddo
    enddo
  enddo

END SUBROUTINE swap_kij_to_ikj

! ****************************************************************************

  subroutine combinit_user(temp_veg_class, temp_glat, &
              temp_glon, nfuel_cat)

  
  use mem_grid, only: &
       jdim,          & ! intent(in)
       ngrid,         & ! intent(in)
       polelat,       & ! intent(in)
       polelon,       & ! intent(in)
       centlat,       & ! intent(in)
       centlon,       & ! intent(in)
       deltax,        & ! intent(in)
       deltay,        & ! intent(in)
       nstratx,       & ! intent(in)
       nstraty,       & ! intent(in)
       nxtnest,       & ! intent(in)
       nnxp,          & ! intent(in)
       nnyp,          & ! intent(in)
       ninest,        & ! intent(inout)
       njnest,        & ! intent(inout)
       deltaxn,       & ! intent(out)
       deltayn,       & ! intent(out)
       platn,         & ! intent(out)
       plonn,         & ! intent(out)
    !   xmn,          & ! intet(in)
    !   ymn,            & !intente(in)
       npatch
   !    xtn,           & ! intent(in)
  !     ytn             ! intent(in)
    
 ! use rconstants, only: erad        
  use module_model_constants, only: DEGRAD
  
  implicit none

  
  integer, parameter:: escrever = 1 !FLAG para escrever no ficheiro para ser lido no grads (1 => sim)
  
  real, intent(OUT) :: nfuel_cat (nnxp(1),nnyp(1))

  integer :: b,i,j  
  !integer, parameter :: nlinhas = 24286 !YDEF
  !integer,parameter :: ncolun = 16439 !XDEF
  !real,parameter :: distan = 0.00022511079999887
  !real,parameter :: lat0 = 36.855458166893
  !real,parameter :: lon0 = -9.7545209749261
  character(len=255) :: filename
  logical :: file_exists
  logical :: value_unit
  integer :: iosize
  character(len=80)::fname
  real, INTENT(IN) ::temp_veg_class(nnxp(1),nnyp(1),npatch)
  real, dimension(nnxp(1),nnyp(1)) :: temp_glat, temp_glon
  character(len=80) :: dummy_str

     !print*,"sfire -estou na rotina dos combustiveis"

     !MFF- Build filename
     write(dummy_str,"(F8.1)")  deltay
     filename = "nfuel_cat_" // adjustl(dummy_str)

     write(dummy_str,"(F8.1)")  deltax
     filename = trim(filename) // "_" // adjustl(dummy_str)

     write(dummy_str,"(a,i0,a,i0,a4)") "_", nnyp(1),"_",nnxp(1),".bin"
     filename = trim(filename) // trim(adjustl(dummy_str)) 
        
     inquire(file=filename, exist=file_exists)
     if(file_exists) then
          print*,"Lendo o ficheiro ",filename
          b=22
          inquire(b, opened=value_unit)
          print*, value_unit
          
          open(b, file=trim(adjustl(filename)), status='old', &
               access='stream', form='unformatted', action='read')
          
          READ(b) ((nfuel_cat(i,j),i=1,nnxp(1)),j=1,nnyp(1))
          close(b)
 
     else
          print*,"Nao existe o ficheiro ",filename
          print*,"O programa escreverá glat_glon.bin e sairá"

          b=22
          inquire(b, opened=value_unit)
          print*, value_unit
          write(fname,3)'glat_glon.bin'
          3 format(a18)
          open(b,file=trim(adjustl(fname)), form='unformatted', action='write', &
                    status='replace', access='stream')
                    
          WRITE(b) nnyp(1),nnxp(1),deltay,deltax
          WRITE(b) (((temp_glat(i,j)),i=1,nnxp(1)),j=1,nnyp(1))
          WRITE(b) (((temp_glon(i,j)),i=1,nnxp(1)),j=1,nnyp(1))
          
          close(b)
 !         call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
           stop -1
 
    endif

   print*,"sfire -vou sair da rotina dos combustiveis"   
   return
  
  end subroutine combinit_user
  
  SUBROUTINE get_ijk_from_subgrid ( config_flags )
    USE mem_grid ,  ONLY : nnzp,nnxp,nnyp 
    USE ModNamelistsfireFile 
    
    TYPE (grid_config_rec_type), TARGET :: config_flags 
    
   

 !- converting WRF setting to BRAMS
     !alterei o indice inicial ids e jds
    config_flags%ids=3   ;config_flags%ide=nnxp(1)-2 ;config_flags%jds=3   ;config_flags%jde=nnyp(1)-2 ;config_flags%kds=2; config_flags%kde=nnzp(1)	      
      config_flags%ims=1   ;config_flags%ime=nnxp(1) ;config_flags%jms=1   ;config_flags%jme=nnyp(1) ;config_flags%kms=1; config_flags%kme=nnzp(1)
     ! config_flags%ids=3   ;config_flags%ide=mxp-2 ;config_flags%jds=3   ;config_flags%jde=myp-2 ;config_flags%kds=1; config_flags%kde=mzp	      
    !  config_flags%ims=1   ;config_flags%ime=mxp ;config_flags%jms=1   ;config_flags%jme=myp ;config_flags%kms=1; config_flags%kme=mzp			
      
      !config_flags%ips=ia+1;config_flags%ipe=iz-2;config_flags%jps=ja+1;config_flags%jpe=jz-2;config_flags%kps=1; config_flags%kpe=mzp			  
      !config_flags%its=ia  ;config_flags%ite=iz  ;config_flags%jts=ja  ;config_flags%jte=jz  ;config_flags%kts=1; config_flags%kte=mzp-1
      
     
        ! substitui pela de baixo:
	
       config_flags%ips=config_flags%ids;config_flags%ipe=config_flags%ide;config_flags%jps=config_flags%jds;
	config_flags%jpe=config_flags%jde;config_flags%kps=2;config_flags%kpe=nnzp(1)			  
     !config_flags%ips=config_flags%ids;config_flags%ipe=config_flags%ide;config_flags%jps=config_flags%jds;
	!config_flags%jpe=config_flags%jde;config_flags%kps=1;config_flags%kpe=mzp
    
   print*,"PROG::(ids,ide,jds,jde,kds,kde)=",config_flags%ids,config_flags%ide,&
   config_flags%jds,config_flags%jde,config_flags%kds,config_flags%kde
   print*,"PROG::(ims,ime,jms,jme,kms,kme)=",config_flags%ims,config_flags%ime,&
   config_flags%jms,config_flags%jme,config_flags%kms,config_flags%kme
   print*,"PROG::(ips,ipe,jps,jpe,kps,kpe)=",config_flags%ips,config_flags%ipe,&
   config_flags%jps,config_flags%jpe,config_flags%kps,config_flags%kpe 
     
     config_flags%ifds = config_flags%ids
     config_flags%ifde = config_flags%ide * config_flags%sr_x
     config_flags%jfds = config_flags%jds
     config_flags%jfde = config_flags%jde * config_flags%sr_y
     config_flags%kfds = config_flags%kds 
     config_flags%kfde = config_flags%kde
                        
     config_flags%ifms = (config_flags%ims-1)*config_flags%sr_x+1
     config_flags%ifme = config_flags%ime * config_flags%sr_x
     config_flags%jfms = (config_flags%jms-1)*config_flags%sr_y+1
     config_flags%jfme = config_flags%jme * config_flags%sr_y
     config_flags%kfms = config_flags%kms
     config_flags%kfme = config_flags%kme 
                    
     config_flags%ifps = (config_flags%ips-3)*config_flags%sr_x+3
     !config_flags%ifps = (config_flags%ips-1)*config_flags%sr_x+1
     config_flags%ifpe = config_flags%ipe * config_flags%sr_x
     config_flags%jfps = (config_flags%jps-3)*config_flags%sr_y+3
     !config_flags%jfps = (config_flags%jps-1)*config_flags%sr_y+1
     config_flags%jfpe =  config_flags%jpe * config_flags%sr_y
     config_flags%kfps = config_flags%kps
     config_flags%kfpe = config_flags%kpe
     
   RETURN
   END SUBROUTINE get_ijk_from_subgrid



   subroutine GlobalLatLonGrid(temp_glon,temp_glat)

    use mem_grid, only: nnxp, nnyp,ngrid 
                
   
    
    implicit none

    integer :: i
    integer :: ierr
    real, INTENT(IN) ::temp_glat(nnxp(ngrid),nnyp(ngrid))
    real, INTENT(IN) ::temp_glon(nnxp(ngrid),nnyp(ngrid))
    real :: delLon, delLat
    real :: deltaSum
    integer :: nLon, nLat
    real :: firstLon, firstLat
    real :: latMin, latMax, latMinBrams, latMaxBrams
    real :: lonMin, lonMax, lonMinBrams, lonMaxBrams
    real :: loni(ngrid)
    real :: lonf(ngrid)
    real :: lati(ngrid)
    real :: latf(ngrid)
    character(len=8) :: c0
    character(len=*), parameter :: h="**(GlobalLatLonGrid)**"
    logical, parameter :: dumpLocal=.true. 
    logical, parameter :: project = .true.
    real, allocatable :: lon(:)    ! longitudes
    real, allocatable :: lat(:)    ! longitudes
    
    
    print*, "DEBUG GlobalLatLonGrid :: ", h

    deltaSum = sum&
         (&
         (temp_glon(nnxp(ngrid),:) - temp_glon(1,:))&
         /(nnxp(ngrid)-1)&
         )
    delLon = deltaSum / nnyp(ngrid)
    print*, "DEBUG GlobalLatLonGrid - delLon ::", delLon
    ! same procedure for latitudes:
    ! sum of average delta latitude over all longitudes
    ! divided by number of longitudes

    deltaSum = sum&
         (&
         (temp_glat(:,nnyp(ngrid)) - temp_glat(:,1))&
         /(nnyp(ngrid)-1)&
         )
    delLat = deltaSum / nnxp(ngrid)
    print*, "DEBUG GlobalLatLonGrid - delLat ::", delLat
    ! Step 2: Compute origin and number of points, restricted to
    ! user selected area and projection  specified at namelist file

    if (project) then

       ! if projection:
       !   find BRAMS grid extremes
       !   intersect with user selected area
       !   define origin and number of points keeping original delta

       ! envelope BRAMS grid area at earth's surface

       lonMinBrams = minval(temp_glon)
       lonMaxBrams = maxval(temp_glon)
       latMinBrams = minval(temp_glat)
       latMaxBrams = maxval(temp_glat)
       
       print*, "DEBUG:: lonMinBrams = minval(temp_glon)",lonMinBrams
       print*, "DEBUG:: lonMaxBrams = maxval(temp_glon)",lonMaxBrams
       print*, "DEBUG:: latMinBrams = minval(temp_glat)",latMinBrams
       print*, "DEBUG:: latMaxBrams = maxval(temp_glat)",latMaxBrams
       
       ! intersection of BRAMS envelope and user defined area
       loni(ngrid) = -180.
       lonf(ngrid) = 180.
       lati(ngrid) = -90.
       latf(ngrid) = 90.
       
       !ATENCAO::: VER SE MUDAMOS MAIS TARDE PARA ISTO::oneNamelistFile%loni(ngrid)), oneNamelistFile%lonf(ngrid)),oneNamelistFile%lati(ngrid)), oneNamelistFile%latf(ngrid))
       lonMin = max(lonMinBrams, loni(ngrid))
       lonMax = min(lonMaxBrams, lonf(ngrid))
       latMin = max(latMinBrams, lati(ngrid))
       latMax = min(latMaxBrams, latf(ngrid))
       print*, "DEBUG:: lonMin =", lonMin
       print*, "DEBUG:: lonMax =", lonMax
       print*, "DEBUG:: latMin =", latMin
       print*, "DEBUG:: latMax =", latMax
       ! post grid origin

       firstLon = lonMin
       firstLat = latMin

       ! post grid number of points (#intervals + 1)

       nLon = 1 + &
            ceiling((lonMax-lonMin)/delLon)
       nLat = 1 + &
            ceiling((latMax-latMin)/delLat)
     
     print*,"nlon,nlat - project",nLon,nLat,project
    else

       ! if no projection:
       !   find BRAMS grid origin as average first latitude and average first longitude
       !   compute last latitude and last longitude keeping delta and # points
       !   intersect with user selected area
       !   define origin and number of points keeping original delta

       ! intersection of BRAMS envelope and user defined area
       loni(ngrid) = -180.
       lonf(ngrid) = 180.
       lati(ngrid) = -90.
       latf(ngrid) = 90.

       lonMinBrams = sum(temp_glon(1,:))/nnyp(ngrid)
       latMinBrams = sum(temp_glat(:,1))/nnxp(ngrid)

       lonMaxBrams = lonMinBrams + delLon*real(nnxp(ngrid)-1)
       latMaxBrams = latMinBrams + delLat*real(nnyp(ngrid)-1)


       !ATENCAO::: VER SE MUDAMOS MAIS TARDE PARA ISTO::oneNamelistFile%loni(ngrid)), oneNamelistFile%lonf(ngrid)),oneNamelistFile%lati(ngrid)), oneNamelistFile%latf(ngrid))
       lonMin = max(lonMinBrams, loni(ngrid))
       lonMax = min(lonMaxBrams, lonf(ngrid))
       latMin = max(latMinBrams, lati(ngrid))
       latMax = min(latMaxBrams, latf(ngrid))
       firstLon = lonMin
       firstLat = latMin
       print*,"lonMin, lonMax,latMin,latMax =",lonMin, lonMax,latMin,latMax
       print*,"firstLon, firstLat =",firstLon, firstLat
       if ( (lonMin == lonMinBrams) .and. &
            (lonMax == lonMaxBrams) .and. &
            (latMin == latMinBrams) .and. &
            (latMax == latMaxBrams)         ) then

          nLon = nnxp(ngrid)
          nLat = nnyp(ngrid)
	  print*,"vim pelo NNXP E NNYP coorde"
       else
          print*,"VOU PELO CEILING"
          nLon = 1 + &
               ceiling((lonMax-lonMin)/delLon)
          nLat = 1 + &
               ceiling((latMax-latMin)/delLat)
          print*,"nLon ,nLat =",nLon ,nLat
       end if
    end if

    ! Step 3: Compute longitude and latitude points
    ! given first point, delta and number of points

    allocate(lon(nLon), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") nLon
       call fatal_error(h//" allocate lon("//&
            trim(adjustl(c0))//") fails")
    end if

    call Axis(nLon, firstLon, &
         delLon, lon)

    allocate(lat(nLat), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") nLat
       call fatal_error(h//" allocate lat("//&
            trim(adjustl(c0))//") fails")
    end if 
    

    call Axis(nLat, firstLat, &
         delLat, lat)
 

    if (dumpLocal) then
       open(15,file='coordenadas.dat',form='formatted',status='replace')
       write(15,"(a6,i8)") "nLon =",nLon
       write(15,"(a6,i8)") "nLat =", nLat
       write(15,"(a10,f12.7)") "firstLon =", firstLon
       write(15,"(f12.7)") lon
       write(15,"(a8,f12.7)") "delLon =",delLon
       write(15,"(a10,f12.7)") "firstLat =",firstLat 
       write(15,"(f12.7)") lat
       write(15,"(a8,f12.7)") "delLat =", delLat
       close(15)
       !call MsgDump (h//" lat-lon has size"//&
       !     " ("//trim(adjustl(c0))//","//trim(adjustl(c1))//")"//&
       !     " with lon,lat = "//&
       !     " ("//trim(adjustl(c2))//":"//trim(adjustl(c3))//":"//trim(adjustl(c4))//", "//&
       !     trim(adjustl(c5))//":"//trim(adjustl(c6))//":"//trim(adjustl(c7))//")")
       
    end if
  end subroutine GlobalLatLonGrid

  subroutine Axis(nVal, firstVal, delVal, val)
    integer, intent(in) :: nVal
    real, intent(in) :: firstVal
    real, intent(in) :: delVal
    real, intent(out) :: val(nVal)

    integer :: i

    print*, 'DEBUG :: func axis'
    do i = 1, nVal
       val(i) = firstVal + real(i-1)*delVal
    end do
    
    return
  end subroutine Axis



  
  
