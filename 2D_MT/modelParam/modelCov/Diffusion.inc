!-------------------------------------------------------------------!
! 2D_MT model covariance by Weerachai Siripunvaraporn: procedures.  !
!-------------------------------------------------------------------!
  function multBy_Cm(m1) result (m2)

   ! Multiplies by the square root of the model covariance,
   ! which is viewed as a smoothing operator. Intended
   ! to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
   ! For efficiency, CmSqrt is a saved, private variable inside
   ! the modelParam module. Before this routine can be called,
   ! it has to be initialized by calling create_CmSqrt(m).

    type (modelParam_t), intent(in)     :: m1
    type (modelParam_t)				    :: m2

    if (.not. CmSqrt%allocated) then
    	call create_CmSqrt(m1)
    end if

    m2 = m1
    call solveDiff(CmSqrt,m2)
    m2%temporary = .true.

  end function multBy_Cm
! *******************************************************************
  function multBy_CmSqrt(m1) result (m2)

   ! Multiplies by the square root of the model covariance,
   ! which is viewed as a smoothing operator. Intended
   ! to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
   ! For efficiency, CmSqrt is a saved, private variable inside
   ! the modelParam module. Before this routine can be called,
   ! it has to be initialized by calling create_CmSqrt(m).

    type (modelParam_t), intent(in)     :: m1
    type (modelParam_t)				    :: m2

    if (.not. CmSqrt%allocated) then
    	call create_CmSqrt(m1)
    end if

    m2 = m1
    call solveDiff(CmSqrt,m2)
    m2%temporary = .true.

  end function multBy_CmSqrt

! *******************************************************************
  subroutine create_CmSqrt(m,cfile)

    type (modelParam_t), intent(in)     :: m
    character(*), intent(in), optional  :: cfile

    if (present(cfile)) then
      ! do nothing: a placeholder for other covariance operators
    end if

    call setup1DCM(m%grid,CmSqrt)

    CmSqrt%allocated = .true.

  end subroutine create_CmSqrt

! *******************************************************************
  subroutine deall_CmSqrt()

    integer                             :: istat

    deallocate(CmSqrt%YDIF, STAT = istat)
    deallocate(CmSqrt%ZDIF, STAT = istat)

  end subroutine deall_CmSqrt

! *******************************************************************
!
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.
!
!----------------------------------------------------------------------
!
!     SOLVING DIFFUSION EQU. IN EXPLICIT WAY.
!...  Copyright (C) 2006 Weerachai Siripunvaraporn. All rights reserved.


  subroutine setup1DCM(gd,Cm)

  type (grid_t),   intent(in)    :: gd
  type (modelCov_t), intent(out)   :: Cm

  real (kind=prec), pointer, dimension(:,:)   ::  YDF
  real (kind=prec), pointer, dimension(:,:)   ::  ZDF

  integer  :: iy, iz, NzEarth ,Ny, istat
  real (kind=prec) :: ylen, zlen

  Ny      = gd%Ny
  NzEarth = gd%Nz - gd%Nza
  allocate(Cm%YDIF(2,gd%Ny,NzEarth),STAT = istat)
  allocate(Cm%ZDIF(2,gd%Ny,NzEarth),STAT = istat)
  Cm%ylen = 0.05
  Cm%zlen = 0.05

  ylen = Cm%ylen
  zlen = Cm%zlen

  ! horizontal diffusion
  do iz = 1,NzEarth
    call setupYCM(Ny,ylen,YDF)
    do iy = 1,Ny
      Cm%YDIF(1,iy,iz) = YDF(1,iy)
      Cm%YDIF(2,iy,iz) = YDF(2,iy)
    enddo
  	deallocate(YDF,STAT=istat)
  enddo

  ! vertical diffusion
  do iy = 1,Ny
    call setupZCM(NzEarth,zlen,ZDF)
    do iz = 1,NzEarth
      Cm%ZDIF(1,iy,iz) = ZDF(1,iz)
      Cm%ZDIF(2,iy,iz) = ZDF(2,iz)
    enddo
    deallocate(ZDF,STAT=istat)
  enddo

  end subroutine setup1DCM

!------------------------------------------------------------------------------

  subroutine setupYCM(Ny,ylen,YDF)

  !use wsLAPACK
  real (kind=prec), pointer, dimension(:,:)   ::  YDF
  integer  :: Ny
  real (kind=prec) :: ylen

  integer  :: iy,info,istat
  real (kind=prec) :: cr, cl
  real (kind=prec), allocatable, dimension(:)   ::  dfyc, dfyr


  allocate(YDF(2,Ny),STAT=istat)
  allocate(dfyc(Ny))
  allocate(dfyr(Ny))
  do iy = 1,Ny
    dfyc(iy) = 0.
    dfyr(iy) = 0.
  enddo

  do iy = 1,Ny
    cr = ylen
    cl = ylen

    dfyc(iy) = 1. + cr + cl
    dfyr(iy) = -cr
!   dfyc(iy) = 1. - cr - cl
!   dfyr(iy) = cr
  enddo

  do iy =1,Ny
    YDF(1,iy) = 0.
    YDF(2,iy) = 0.
  enddo

  do iy =2,Ny
    YDF(1,iy) = dfyr(iy-1)
  enddo
  do iy =1,Ny
    YDF(2,iy) = dfyc(iy)
  enddo

  call dpbtrf('U',Ny,1,YDF,2,info)

  deallocate(dfyc)
  deallocate(dfyr)

  end subroutine setupYCM

!------------------------------------------------------------------------------

  subroutine setupZCM(Nzb,zlen,ZDF)

  !use wsLAPACK
  real (kind=prec), pointer,  dimension(:,:)   ::  ZDF
  integer  :: Nzb
  real (kind=prec) :: zlen

  integer  :: iz,info,istat
  real (kind=prec) :: cr, cl
  real (kind=prec), allocatable, dimension(:)   ::  dfzc, dfzr


  allocate(ZDF(2,Nzb),STAT=istat)
  allocate(dfzc(Nzb))
  allocate(dfzr(Nzb))
  do iz = 1,Nzb
    dfzc(iz) = 0.
    dfzr(iz) = 0.
  enddo

  do iz = 1,Nzb
    cr = zlen
    cl = zlen

    dfzc(iz) = 1. + cr + cl
    dfzr(iz) = -cr
!   dfzc(iz) = 1. - cr - cl
!   dfzr(iz) =  cr
  enddo

  do iz =1,Nzb
    ZDF(1,iz) = 0.
    ZDF(2,iz) = 0.
  enddo

  do iz =2,Nzb
    ZDF(1,iz) = dfzr(iz-1)
  enddo
  do iz =1,Nzb
    ZDF(2,iz) = dfzc(iz)
  enddo

  call dpbtrf('U',Nzb,1,ZDF,2,info)

  deallocate(dfzc)
  deallocate(dfzr)

  end subroutine setupZCM

!------------------------------------------------------------------------------

  subroutine solveDiff(Cm,CmJ)

  type (modelCov_t),   intent(in)    :: Cm
  type (modelParam_t), intent(inout) :: CmJ

  integer            :: Ny,Nz,it, dTime

  Ny = CmJ%Ny
  Nz = CmJ%NzEarth

  dTime = 10

! do iz = 1,Nz
!   do iy = 1,Ny
!     write(45,*) CmJ%v(iy,iz)
!   enddo
! enddo

  do it = 1,dTime
    call yDiff(Cm,CmJ)
    call zDiff(Cm,CmJ)
  enddo
  do it = 1,dTime
    call zDiff(Cm,CmJ)
    call yDiff(Cm,CmJ)
  enddo

! do iz = 1,Nz
!   do iy = 1,Ny
!     write(46,*) CmJ%v(iy,iz)
!   enddo
! enddo

  end subroutine solveDiff

!------------------------------------------------------------------------------

  subroutine yDiff(Cm,CmJ)

  !use wsLAPACK
  type (modelCov_t),   intent(in)    :: Cm
  type (modelParam_t), intent(inout) :: CmJ

  real (kind=prec), allocatable, dimension(:,:)   ::  YDF
  real (kind=prec), allocatable, dimension(:)     ::  uy

  integer            :: Ny,Nz,iy,iz,info,istat

  Ny = CmJ%Ny
  Nz = CmJ%NzEarth

  allocate(YDF(2,Ny))
  allocate(uy(Ny))

  do iz = 1,Nz
    do iy = 1,Ny
      uy(iy)    = CmJ%v(iy,iz)
      YDF(1,iy) = Cm%YDIF(1,iy,iz)
      YDF(2,iy) = Cm%YDIF(2,iy,iz)
    enddo

    call dpbtrs('U',Ny,1,1,YDF,2,uy,Ny,info)

    do iy = 1,Ny
      CmJ%v(iy,iz) = uy(iy)
    enddo

  enddo

  deallocate(YDF,STAT=istat)
  deallocate(uy,STAT=istat)

  end subroutine yDiff

!------------------------------------------------------------------------------

  subroutine zDiff(Cm,CmJ)

  !use wsLAPACK
  type (modelCov_t),   intent(in)    :: Cm
  type (modelParam_t), intent(inout) :: CmJ

  real (kind=prec), allocatable, dimension(:,:)   ::  ZDF
  real (kind=prec), allocatable, dimension(:)     ::  uz


  integer            :: Ny,Nz,iy,iz,info,istat

  Ny = CmJ%Ny
  Nz = CmJ%NzEarth

  allocate(ZDF(2,Nz))
  allocate(uz(Nz))

  do iy = 1,Ny
    do iz = 1,Nz
      uz(iz)    = CmJ%v(iy,iz)
      ZDF(1,iz) = Cm%ZDIF(1,iy,iz)
      ZDF(2,iz) = Cm%ZDIF(2,iy,iz)
    enddo

    call dpbtrs('U',Nz,1,1,ZDF,2,uz,Nz,info)

    do iz = 1,Nz
      CmJ%v(iy,iz) = uz(iz)
    enddo

  enddo

  deallocate(ZDF,STAT=istat)
  deallocate(uz,STAT=istat)

  end subroutine zDiff
