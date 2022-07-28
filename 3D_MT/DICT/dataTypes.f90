! *****************************************************************************
module dataTypes
  ! This module contains the data type dictionary (typeDict) for 3D MT

  use math_constants
  use utilities

  implicit none

  public			:: setup_typeDict, deall_typeDict

  !  stores information about the "data type" -- which could include
  !   information that is relevant to transmitter, receiver, or both
  !   E.g., for 2D MT, "mode" (TE/TM) is relevant to both receiver and
  !    transmitter.  For 3D MT, this can be used to distinguish between
  !    full impedance tensors, impedance tensors+vertical field TFs,
  !    off-diagonal impedance only, interstation TFs, etc.
  !    Different data types may correspond to different EM solutions
  !     (TE vs. TM; joint inversion of data from multiple geophysical
  !        techniques), or not (different receiver configurations in 3D MT).
  !    In some cases, multiple transmitter dictionaries may be needed (and
  !     which to use would be determined in, e.g., EMsolve).
  !    Similarly, in some cases multiple receiver dictionaries may or may not
  !     be required when there are multiple dataTypes
  !       (e.g., interstation TFs require more information, to
  !        specify a pair of site locations; adding more local TF components
  !         would not require multiple receiver dictionaries)
  type :: dataType

     ! The following two attributes (essentially receiver attributes)
     ! must be defined for all data types.  These are accessed
     ! and used by the top-level inversion routines.
     !
     ! isComplex is also stored in the dataVector; the sole purpose of the dataType
     ! should be to identify the mapping from the EM solution to the data vector.
     logical           :: isComplex = .false.

     ! a user-friendly description of this data type
     character(80)     :: name = ''

     ! for 3D MT data types will initially be used only to distinguish
     ! between local transfer function types (later maybe add interstation TFs)
     integer		   :: tfType

     ! the units of the data type. If a value has different units, it's a new data type.
     character(80)     :: units

     ! the number of components in the data type
     integer           :: nComp

     ! the (real or complex) data type components in a fixed order as given;
     ! the file can have a different component order on input.
     character(15), pointer, dimension(:) :: id

     ! these lists contain the indices into the data vector for each data type;
     ! they make it possible to sort the data by receiver for output.
     ! no data denoted by zero index; dimensions (nTx) and (nTx,nRx).
     integer, pointer, dimension(:)   :: tx_index
     integer, pointer, dimension(:)   :: dt_index
     integer, pointer, dimension(:,:) :: rx_index

  end type dataType

  ! data type dictionary must be public; some attributes are referenced
  !   by top-level inversion routines
  type (dataType), pointer, save, public, dimension(:) :: typeDict

  ! For impedance plus vertical components, use data types 1 & 3.
  ! Rho & Phase added on behalf of Kristina Tietze, GFZ-Potsdam by Naser Meqbel
  integer, parameter   :: Full_Impedance = 1
  integer, parameter   :: Off_Diagonal_Impedance = 2
  integer, parameter   :: Full_Vertical_Components = 3
  integer, parameter   :: Full_Interstation_TF = 4
  integer, parameter   :: Off_Diagonal_Rho_Phase = 5
  integer, parameter   :: Phase_Tensor = 6     
  !Add EM_Fileds                                                WangPeijie 2020.12.03
  integer, parameter   :: EM_Fields = 7
 

Contains


!**************************************************************************
! Initializes and sets up data type dictionary
  subroutine setup_typeDict()

  	 integer     :: istat
                 
    ! allocate(typeDict(6),STAT=istat)       
    !typeDict(6)->typeDict(7) EM_Fileds             WangPeijie 2020.12.03                
     allocate(typeDict(7),STAT=istat)

     typeDict(Full_Impedance)%name = 'Full_Impedance'
     typeDict(Full_Impedance)%isComplex = .true.
     typeDict(Full_Impedance)%tfType    = Full_Impedance
     typeDict(Full_Impedance)%units     = '[V/m]/[T]'
     typeDict(Full_Impedance)%nComp     = 8
     allocate(typeDict(Full_Impedance)%id(4),STAT=istat)
     typeDict(Full_Impedance)%id(1)    = 'ZXX'
     typeDict(Full_Impedance)%id(2)    = 'ZXY'
     typeDict(Full_Impedance)%id(3)    = 'ZYX'
     typeDict(Full_Impedance)%id(4)    = 'ZYY'

     typeDict(Off_Diagonal_Impedance)%name = 'Off_Diagonal_Impedance'
     typeDict(Off_Diagonal_Impedance)%isComplex = .true.
     typeDict(Off_Diagonal_Impedance)%tfType    = Off_Diagonal_Impedance
     typeDict(Off_Diagonal_Impedance)%units     = '[V/m]/[T]'
     typeDict(Off_Diagonal_Impedance)%nComp     = 4
     allocate(typeDict(Off_Diagonal_Impedance)%id(2),STAT=istat)
     typeDict(Off_Diagonal_Impedance)%id(1)    = 'ZXY'
     typeDict(Off_Diagonal_Impedance)%id(2)    = 'ZYX'

     typeDict(Full_Vertical_Components)%name = 'Full_Vertical_Components'
     typeDict(Full_Vertical_Components)%isComplex = .true.
     typeDict(Full_Vertical_Components)%tfType    = Full_Vertical_Components
     typeDict(Full_Vertical_Components)%units     = '[]'
     typeDict(Full_Vertical_Components)%nComp     = 4
     allocate(typeDict(Full_Vertical_Components)%id(2),STAT=istat)
     typeDict(Full_Vertical_Components)%id(1)    = 'TX '
     typeDict(Full_Vertical_Components)%id(2)    = 'TY '

     typeDict(Full_Interstation_TF)%name = 'Full_Interstation_TF'
     typeDict(Full_Interstation_TF)%isComplex = .true.
     typeDict(Full_Interstation_TF)%tfType    = Full_Interstation_TF
     typeDict(Full_Interstation_TF)%units     = '[]'
     typeDict(Full_Interstation_TF)%nComp     = 8
     allocate(typeDict(Full_Interstation_TF)%id(4),STAT=istat)
     typeDict(Full_Interstation_TF)%id(1)    = 'MXX'
     typeDict(Full_Interstation_TF)%id(2)    = 'MXY'
     typeDict(Full_Interstation_TF)%id(3)    = 'MYX'
     typeDict(Full_Interstation_TF)%id(4)    = 'MYY'

 	 typeDict(Off_Diagonal_Rho_Phase)%name = 'Off_Diagonal_Rho_Phase'
     typeDict(Off_Diagonal_Rho_Phase)%isComplex = .false.
     typeDict(Off_Diagonal_Rho_Phase)%tfType    = Off_Diagonal_Rho_Phase
     typeDict(Off_Diagonal_Rho_Phase)%units     = '[]'
     typeDict(Off_Diagonal_Rho_Phase)%nComp     = 4
     allocate(typeDict(Off_Diagonal_Rho_Phase)%id(4),STAT=istat)
     typeDict(Off_Diagonal_Rho_Phase)%id(1) = 'RHOXY'
     typeDict(Off_Diagonal_Rho_Phase)%id(2) = 'PHSXY'
     typeDict(Off_Diagonal_Rho_Phase)%id(3) = 'RHOYX'
     typeDict(Off_Diagonal_Rho_Phase)%id(4) = 'PHSYX'
	 	 
	 typeDict(Phase_Tensor)%name = 'Phase_Tensor'
     typeDict(Phase_Tensor)%isComplex = .false.
     typeDict(Phase_Tensor)%tfType    = Phase_Tensor
     typeDict(Phase_Tensor)%units     = '[]'
     typeDict(Phase_Tensor)%nComp     = 4
     allocate(typeDict(Phase_Tensor)%id(4),STAT=istat)
     typeDict(Phase_Tensor)%id(1) = 'PTXX'
     typeDict(Phase_Tensor)%id(2) = 'PTXY'
     typeDict(Phase_Tensor)%id(3) = 'PTYX'
     typeDict(Phase_Tensor)%id(4) = 'PTYY'       
	 	                                                                     
  !Add EM_Fileds                                                WangPeijie 2020.12.03    
	 typeDict(EM_Fields)%name = 'EM_Fields'
     typeDict(EM_Fields)%isComplex = .true.
     typeDict(EM_Fields)%tfType    = EM_Fields
     typeDict(EM_Fields)%units     = '[V/m]/[T]'
     typeDict(EM_Fields)%nComp     = 20
     allocate(typeDict(EM_Fields)%id(10),STAT=istat)
     typeDict(EM_Fields)%id(1) = 'EPolEx'
     typeDict(EM_Fields)%id(2) = 'EPolEy'
     typeDict(EM_Fields)%id(3) = 'EPolHx'
     typeDict(EM_Fields)%id(4) = 'EPolHy'   
     typeDict(EM_Fields)%id(5) = 'EPolHz'       
     typeDict(EM_Fields)%id(6) = 'BPolEx'
     typeDict(EM_Fields)%id(7) = 'BPolEy'
     typeDict(EM_Fields)%id(8) = 'BPolHx'
     typeDict(EM_Fields)%id(9) = 'BPolHy'   
     typeDict(EM_Fields)%id(10) = 'BPolHz'
	 

  end subroutine setup_typeDict

! **************************************************************************
! Cleans up and deletes type dictionary at end of program execution
  subroutine deall_typeDict()

	integer     :: j, istat

	if (associated(typeDict)) then

	   do j = 1,size(typeDict)
	      if (associated(typeDict(j)%id)) then
	         deallocate(typeDict(j)%id,STAT=istat)
	      end if
          if (associated(typeDict(j)%tx_index)) then
             deallocate(typeDict(j)%tx_index,STAT=istat)
          end if
          if (associated(typeDict(j)%dt_index)) then
             deallocate(typeDict(j)%dt_index,STAT=istat)
          end if
          if (associated(typeDict(j)%rx_index)) then
             deallocate(typeDict(j)%rx_index,STAT=istat)
          end if
	   end do

       deallocate(typeDict,STAT=istat)

    end if

  end subroutine deall_typeDict

!**********************************************************************
! Computes the value by which the data must be multiplied to convert
! from the old units to the new units.
! The units may be any of the following.
! 1) SI units for E/B: [V/m]/[T] (used in ModEM code)
! 2) practical units for E/B: [mV/km]/[nT]
! 3) SI units for E/H: [V/m]/[A/m] = Ohm

  function ImpUnits(oldUnits,newUnits) result (SI_factor)
	character(*), intent(in)    :: oldUnits, newUnits
	real(kind=prec)             :: SI_factor
	! local
	real(kind=prec)             :: factor1, factor2

	! if the quantity is dimensionless, do nothing
	if ((index(oldUnits,'[]')>0) .or. (index(newUnits,'[]')>0)) then
	   SI_factor = ONE
	   return
	end if

		! first convert the old units to [V/m]/[T]
		if (index(oldUnits,'[V/m]/[T]')>0) then
		   ! SI units for E/B
		   factor1 = ONE
		else if (index(oldUnits,'[mV/km]/[nT]')>0) then
		   ! practical units for E/B
		   factor1 = ONE * 1000.0
		else if ((index(oldUnits,'[V/m]/[A/m]')>0) .or. (index(oldUnits,'Ohm')>0)) then
		   ! SI units for E/H
		   factor1 = ONE * 1000.0 * 10000.0/(4*PI) ! approx. 796000.0
		else
		   call errStop('Unknown input units in ImpUnits: '//trim(oldUnits))
		end if

		! now convert [V/m]/[T] to the new units
		if (index(newUnits,'[V/m]/[T]')>0) then
		   ! SI units for E/B
		   factor2 = ONE
		else if (index(newUnits,'[mV/km]/[nT]')>0) then
		   ! practical units for E/B
		   factor2 = ONE / (1000.0)
		else if ((index(newUnits,'[V/m]/[A/m]')>0) .or. (index(newUnits,'Ohm')>0)) then
		   ! SI units for E/H
		   factor2 = ONE / (1000.0 * 10000.0/(4*PI))
		else
		   call errStop('Unknown output units in ImpUnits: '//trim(newUnits))
		end if

		SI_factor = factor1 * factor2


  end function ImpUnits

!**********************************************************************
! Figures out the data type from its name

  function ImpType(typeName) result (dataType)

    character(*), intent(in)    :: typeName
	integer	             	 	:: dataType

    select case (trim(typeName))

       case('Full_Impedance')
          dataType = Full_Impedance

       case('Off_Diagonal_Impedance')
          dataType = Off_Diagonal_Impedance

       case('Full_Vertical_Components')
          dataType = Full_Vertical_Components

       case('Full_Interstation_TF')
          dataType = Full_Interstation_TF

       case('Off_Diagonal_Rho_Phase')
          dataType = Off_Diagonal_Rho_Phase
		  		  
      case('Phase_Tensor')
          dataType = Phase_Tensor     
		  		  
      case('EM_Fields')
          dataType = EM_Fields
		  
       case default
          call errStop('Unknown data type:'//trim(typeName))

    end select

  end function ImpType

!**********************************************************************
! Sorts out the data type header

  function ImpHeader(dataType) result (header)

    integer, intent(in)         :: dataType
    character(200)              :: header

    select case (dataType)

       case(Full_Impedance,Off_Diagonal_Impedance,Full_Vertical_Components)
          header = '# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error'

       case(Full_Interstation_TF)
          header = '# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Ref_Code Ref_Lat Ref_Lon Ref_X(m) Ref_Y(m) Ref_Z(m) Component Real Imag Error'

       case(Off_Diagonal_Rho_Phase,Phase_Tensor)
          header = '# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Value Error'

    end select

  end function ImpHeader

!**********************************************************************
! Figures out the component index from its name for any data type

  function ImpComp(compid,dataType) result (icomp)

    character(*), intent(in)    :: compid
    integer,      intent(in)    :: dataType
    integer                     :: icomp
    ! local
    integer                     :: i, ncomp

    ncomp = typeDict(dataType)%ncomp
    if (typeDict(dataType)%isComplex) then
        ncomp = ncomp/2
    end if
    icomp = 0

    do i = 1,ncomp
        if (index(typeDict(dataType)%id(i),trim(compid))>0) then
            icomp = i
            exit
        end if
    end do

    if (icomp == 0) then
        call errStop('Problem locating the component '//trim(compid)//' in data type '//trim(typeDict(dataType)%name))
    end if

  end function ImpComp

!*************************************************************************
! Supports the old data format file: figures out the data type from header

  function ImpTypeFromHeader(nComp,header) result (dataType)

    integer, intent(in)         :: nComp
    character(*), intent(in)    :: header
    integer                     :: dataType
    ! local
    character(15), allocatable  :: compids(:)
    integer                     :: j,istat

    allocate(compids(nComp),STAT=istat)

    read(header,*) compids

    select case (nComp)
       case(8)
          if (index(compids(1),'Mxx')>0) then
             dataType =  Full_Interstation_TF
          else
             dataType =  Full_Impedance
          end if
       case(4)
          if (index(compids(1),'Tx')>0) then
             dataType =  Full_Vertical_Components
          elseif (index(compids(1),'Zxy')>0) then
             dataType =  Off_Diagonal_Impedance
          elseif (index(compids(1),'Rhoxy')>0) then
             dataType =  Off_Diagonal_Rho_Phase
          end if
    end select

    do j = 1,nComp
        if (compids(j) .ne. typeDict(dataType)%id(j)) then
            call errStop('Wrong order of impedance components in data header')
        end if
    end do

    deallocate(compids,STAT=istat)

  end function ImpTypeFromHeader


!*************************************************************************
! Supports the old data format file: checks the order of components

 subroutine check_header_order(nComp,dataType,header)
 	 integer, intent(in)   	 	 :: nComp
 	 integer, intent(in)   	 	 :: dataType
     character(*), intent(in)    :: header
 	 !Local
 	 character(15), allocatable  :: compids(:)
     integer                     :: j,istat

     allocate(compids(nComp),STAT=istat)
     read(header,*) compids


      do j = 1,nComp
    	if (compids(j) .ne. typeDict(dataType)%id(j)) then
    		call errStop('Wrong order of impedance components in data header')
    	end if
    end do

    deallocate(compids,STAT=istat)

 end subroutine check_header_order


end module dataTypes
