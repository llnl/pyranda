module LES_compact
  USE iso_c_binding
  USE MPI
  !USE LES_input, ONLY : directcom,zerodx,zerody,zerodz,gpu_kernel
  USE LES_stencils
  USE LES_patch, ONLY : patch_type
  USE LES_comm,  ONLY : comm_type,master
  use LES_compact_basetype, only : compact_op1,mesh1_type,comm1_type, &
  																 vbr1x,vbr2x,vbs1x,vbs2x,vbr1y,vbr2y,vbs1y,vbs2y, &
  																 vbr1z,vbr2z,vbs1z,vbs2z,dvopx,dvox,dvopy,dvoy, &
  																 dvopz,dvoz
  use LES_compact_d1, only : compact_op1_d1
  use LES_compact_r3, only : compact_op1_r3
  use LES_compact_r4, only : compact_op1_r4

  IMPLICIT NONE


  logical :: zerodx = .false., zerody = .false., zerodz = .false.
  integer :: gpu_kernel = 1
  integer :: directcom = 1
  
  REAL(KIND=c_double), PARAMETER :: zero=0.0_c_double, one=1.0_c_double
  LOGICAL(c_bool) :: debug=.false.

  TYPE control_type  ! from input or samrai
    LOGICAL(c_bool) :: null_opx=.false.,null_opy=.false.,null_opz=.false.    ! skip operations
    integer(c_int) :: d1spec=1,d2spec=1,d4spec=1,d8spec=1,gfspec=6,sfspec=2,tfspec=8  ! compact scheme
    integer(c_int) :: islspec=1,isrspec=2,imcfspec=1
    INTEGER(c_int) :: directcom = 1  ! parallel solve
  END TYPE control_type
  
  type(control_type) :: control_data

  type compact_type  ! suite of operations
    type(control_type) :: control  ! controls
    integer(c_int) :: mbc(2,3,2), nop(3) ! (boundary, direction, symmetry), number of operators
    real(KIND=c_double) :: dx,dy,dz  ! nomimal grid spacing for derivatives
    type(compact_op1_d1), dimension(2) :: d1x,d1y,d1z
    type(compact_op1_r3), dimension(2) :: d2x,d2y,d2z
    type(compact_op1_r3), dimension(2) :: d4x,d4y,d4z
    type(compact_op1_r4), dimension(2) :: d8x,d8y,d8z
    type(compact_op1_r4), dimension(2) :: gfx,gfy,gfz, sfx,sfy,sfz, tfx,tfy,tfz  ! gaussian, spectral, tophat filters
    type(compact_op1_r3), dimension(2) :: icfx,icfy,icfz
  contains
    procedure :: setup => setup_compact_ops
    procedure :: remove => remove_compact_ops
  end type compact_type
  
  type(compact_type), pointer :: compact_ops => null()

contains

! setup and remove routines

  subroutine setup_compact_ops(cops,patch,com,patch_num,patch_level)
    implicit none
    class(compact_type) :: cops
    type(patch_type), intent(in) :: patch
    type(comm_type), intent(in) :: com
    type(mesh1_type) :: xmsh,ymsh,zmsh
    type(comm1_type) :: xcom,ycom,zcom
    type(compact_weight), pointer :: weight
    logical(c_bool) :: null_op,null_opx,null_opy,null_opz
    logical :: spew = .FALSE.
    integer :: n,spec,patch_num,patch_level !,directcom
    ! get 1D from 3D types in terms of Miranda block and comm types
     xmsh = mesh1_type(patch%ax,patch%nx,patch%dx,patch%x1,patch%xn,patch%bx1,patch%bxn)
     xcom = comm1_type(patch%periodicx,com%xcom,com%xcom_np,com%xcom_id, &
       com%xcom_lo,com%xcom_hi,com%xrange)
     ymsh = mesh1_type(patch%ay,patch%ny,patch%dy,patch%y1,patch%yn,patch%by1,patch%byn)
     ycom = comm1_type(patch%periodicy,com%ycom,com%ycom_np,com%ycom_id, &
       com%ycom_lo,com%ycom_hi,com%yrange)
     zmsh = mesh1_type(patch%az,patch%nz,patch%dz,patch%z1,patch%zn,patch%bz1,patch%bzn)
     zcom = comm1_type(patch%periodicz,com%zcom,com%zcom_np,com%zcom_id, &
       com%zcom_lo,com%zcom_hi,com%zrange)
    ! unknown/periodic (default)
    cops%mbc = 0 ; cops%nop = 1
    ! symmetric-antisymmetric pairs
    if( xmsh%bc1=='SYMM' ) then ; cops%mbc(1,1,1:2) = [1,-1] ; cops%nop(1) = 2 ; endif
    if( xmsh%bcn=='SYMM' ) then ; cops%mbc(2,1,1:2) = [1,-1] ; cops%nop(1) = 2 ; endif
    if( ymsh%bc1=='SYMM' ) then ; cops%mbc(1,2,1:2) = [1,-1] ; cops%nop(2) = 2 ; endif
    if( ymsh%bcn=='SYMM' ) then ; cops%mbc(2,2,1:2) = [1,-1] ; cops%nop(2) = 2 ; endif
    if( zmsh%bc1=='SYMM' ) then ; cops%mbc(1,3,1:2) = [1,-1] ; cops%nop(3) = 2 ; endif
    if( zmsh%bcn=='SYMM' ) then ; cops%mbc(2,3,1:2) = [1,-1] ; cops%nop(3) = 2 ; endif
    ! extended rhs
    if( xmsh%bc1=='XTND' ) cops%mbc(1,1,1:2) = 2
    if( xmsh%bcn=='XTND' ) cops%mbc(2,1,1:2) = 2
    if( ymsh%bc1=='XTND' ) cops%mbc(1,2,1:2) = 2
    if( ymsh%bcn=='XTND' ) cops%mbc(2,2,1:2) = 2
    if( zmsh%bc1=='XTND' ) cops%mbc(1,3,1:2) = 2
    if( zmsh%bcn=='XTND' ) cops%mbc(2,3,1:2) = 2
    ! grid spacing
    cops%dx = xmsh%d ; cops%dy = ymsh%d ; cops%dz = zmsh%d
    ! operation control (from input)
    cops%control%null_opx = ( patch%nx < 4 .or. zerodx )
    cops%control%null_opy = ( patch%ny < 4 .or. zerody )
    cops%control%null_opz = ( patch%nz < 4 .or. zerodz )
    cops%control%directcom = directcom
    ! local copies
    null_opx = cops%control%null_opx ; null_opy = cops%control%null_opy ; null_opz = cops%control%null_opz
    ! make sure weights are set up
!    spew = (com%patcom_id == master) ! .false.
    if( .not. stencils_set ) call setup_stencils(spew)
    ! d1
    spec = cops%control%d1spec
    if( spec < 1 .or. spec > nderiv1 ) then  ! should be a fatal error
      null_op = .true.
      print *,'## invalid d1spec ##'
    else
      null_op = .false.
      weight => compact_weight_d1(spec)
    endif
    do n=1,cops%nop(1)
      cops%d1x(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d1x(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%d1x(n)%directcom = directcom
      !if( spew .and. .not. cops%d1x(n)%null_op ) 
      print *,n,'setting up compact d1x ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%d1y(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d1y(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%d1y(n)%directcom = directcom
      if( spew .and. .not. cops%d1y(n)%null_op ) print *,n,'setting up compact d1y ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%d1z(n)%null_op = null_op
      if( null_op ) cycle
     call cops%d1z(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%d1z(n)%directcom = directcom
      if( spew .and. .not. cops%d1z(n)%null_op ) print *,n,'setting up compact d1z ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! d2
    spec = cops%control%d2spec
    if( spec < 1 .or. spec > nderiv2 ) then
      null_op = .true.
      print *,'## invalid d2spec ##'
    else
      null_op = .false.
      weight => compact_weight_d2(spec)
    endif
    do n=1,cops%nop(1)
      cops%d2x(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d2x(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%d2x(n)%directcom = directcom
      if( spew .and. .not. cops%d2x(n)%null_op ) print *,n,'setting up compact d2x ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%d2y(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d2y(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%d2y(n)%directcom = directcom
      if( spew .and. .not. cops%d2y(n)%null_op ) print *,n,'setting up compact d2y ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%d2z(n)%null_op = null_op
      if( null_op ) cycle
     call cops%d2z(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%d2z(n)%directcom = directcom
      if( spew .and. .not. cops%d2z(n)%null_op ) print *,n,'setting up compact d2z ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! d4
    spec = cops%control%d4spec
    if( spec < 1 .or. spec > nderiv4 ) then
      null_op = .true.
      print *,'## invalid d4spec ##'
    else
      null_op = .false.
      weight => compact_weight_d4(spec)
    endif
    do n=1,cops%nop(1)
      cops%d4x(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d4x(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%d4x(n)%directcom = directcom
      if( spew .and. .not. cops%d4x(n)%null_op ) print *,n,'setting up compact d4x ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%d4y(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d4y(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%d4y(n)%directcom = directcom
      if( spew .and. .not. cops%d4y(n)%null_op ) print *,n,'setting up compact d4y ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%d4z(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d4z(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%d4z(n)%directcom = directcom
      if( spew .and. .not. cops%d4z(n)%null_op ) print *,n,'setting up compact d4z ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! d8
    spec = cops%control%d8spec
    if( spec < 1 .or. spec > nderiv8 ) then
      null_op = .true.
      print *,'## invalid d8spec ##'
    else
      null_op = .false.
      weight => compact_weight_d8(spec)
    endif
    do n=1,cops%nop(1)
      cops%d8x(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d8x(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%d8x(n)%directcom = directcom
      if( spew .and. .not. cops%d8x(n)%null_op ) print *,n,'setting up compact d8x ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%d8y(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d8y(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%d8y(n)%directcom = directcom
      if( spew .and. .not. cops%d8y(n)%null_op ) print *,n,'setting up compact d8y ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%d8z(n)%null_op = null_op
      if( null_op ) cycle
      call cops%d8z(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%d8z(n)%directcom = directcom
      if( spew .and. .not. cops%d8z(n)%null_op ) print *,n,'setting up compact d8z ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! gf
    spec = cops%control%gfspec
    if( spec < 1 .or. spec > nfilter ) then
      null_op = .true.
      print *,'## invalid gfspec ##'
    else
      null_op = .false.
      weight => compact_weight_ff(spec)
    endif
    do n=1,cops%nop(1)
      cops%gfx(n)%null_op = null_op
      if( null_op ) cycle
      call cops%gfx(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%gfx(n)%directcom = directcom
      if( spew .and. .not. cops%gfx(n)%null_op ) print *,n,'setting up compact gfx ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%gfy(n)%null_op = null_op
      if( null_op ) cycle
      call cops%gfy(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%gfy(n)%directcom = directcom
      if( spew .and. .not. cops%gfy(n)%null_op ) print *,n,'setting up compact gfy ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%gfz(n)%null_op = null_op
      if( null_op ) cycle
      call cops%gfz(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%gfz(n)%directcom = directcom
      if( spew .and. .not. cops%gfz(n)%null_op ) print *,n,'setting up compact gfz ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! sf
    spec = cops%control%sfspec
    if( spec < 1 .or. spec > nfilter ) then
      null_op = .true.
      print *,'## invalid sfspec ##'
    else
      null_op = .false.
      weight => compact_weight_ff(spec)
    endif
    do n=1,cops%nop(1)
      cops%sfx(n)%null_op = null_op
      if( null_op ) cycle
      call cops%sfx(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%sfx(n)%directcom = directcom
      if( spew .and. .not. cops%sfx(n)%null_op ) print *,n,'setting up compact sfx ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%sfy(n)%null_op = null_op
      if( null_op ) cycle
      call cops%sfy(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%sfy(n)%directcom = directcom
      if( spew .and. .not. cops%sfy(n)%null_op ) print *,n,'setting up compact sfy ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%sfz(n)%null_op = null_op
      if( null_op ) cycle
      call cops%sfz(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%sfz(n)%directcom = directcom
      if( spew .and. .not. cops%sfz(n)%null_op ) print *,n,'setting up compact sfz ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
    ! tf
    spec = cops%control%tfspec
    if( spec < 1 .or. spec > nfilter ) then
      null_op = .true.
      print *,'## invalid tfspec ##'
    else
      null_op = .false.
      weight => compact_weight_ff(spec)
    endif
    do n=1,cops%nop(1)
      cops%tfx(n)%null_op = null_op
      if( null_op ) cycle
      call cops%tfx(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%tfx(n)%directcom = directcom
      if( spew .and. .not. cops%tfx(n)%null_op ) print *,n,'setting up compact tfx ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%tfy(n)%null_op = null_op
      if( null_op ) cycle
      call cops%tfy(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%tfy(n)%directcom = directcom
      if( spew .and. .not. cops%tfy(n)%null_op ) print *,n,'setting up compact tfy ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%tfz(n)%null_op = null_op
      if( null_op ) cycle
      call cops%tfz(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%tfz(n)%directcom = directcom
      if( spew .and. .not. cops%tfz(n)%null_op ) print *,n,'setting up compact tfz ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do
   
    ! High-order interpolation (zone to face)
    spec = cops%control%imcfspec
    null_op = .false.
    weight => compact_weight_imcf(spec)
    
    do n=1,cops%nop(1)
      cops%icfx(n)%null_op = null_op
      if( null_op ) cycle
      call cops%icfx(n)%setup(weight,xcom,xmsh,cops%mbc(:,1,n),null_opx)
      cops%icfx(n)%directcom = directcom
      !if( spew .and. .not. cops%icfx(n)%null_op ) 
      print *,n,'setting up compact icfx ',weight%description%name,' for ',xmsh%bc1,' ',xmsh%bcn
    end do
    do n=1,cops%nop(2)
      cops%icfy(n)%null_op = null_op
      if( null_op ) cycle
      call cops%icfy(n)%setup(weight,ycom,ymsh,cops%mbc(:,2,n),null_opy)
      cops%icfy(n)%directcom = directcom
      if( spew .and. .not. cops%icfy(n)%null_op ) print *,n,'setting up compact tfy ',weight%description%name,' for ',ymsh%bc1,' ',ymsh%bcn
    end do
    do n=1,cops%nop(3)
      cops%icfz(n)%null_op = null_op
      if( null_op ) cycle
      call cops%icfz(n)%setup(weight,zcom,zmsh,cops%mbc(:,3,n),null_opz)
      cops%icfz(n)%directcom = directcom
      if( spew .and. .not. cops%icfz(n)%null_op ) print *,n,'setting up compact tfz ',weight%description%name,' for ',zmsh%bc1,' ',zmsh%bcn
    end do


    ! Allocate scratchpad for compact operations
    if(patch_level==0 .and. patch_num==0) then	! Am I on the base patch?
      CALL setup_ghost_buffers(patch)
    end if
  end subroutine setup_compact_ops

  subroutine remove_compact_ops(compact,patch_level,patch_num)
    class(compact_type) :: compact
    integer, intent(in) :: patch_num,patch_level
    integer :: n
    ! deallocate buffers on base grid
    if( patch_level==0 .and. patch_num==0 ) CALL remove_ghost_buffers()
    ! free pointer arrays
    do n=1,compact%nop(3)
      if( .not. compact%tfz(n)%null_op ) call compact%tfz(n)%remove()
    end do
    do n=1,compact%nop(2)
      if( .not. compact%tfy(n)%null_op ) call compact%tfy(n)%remove()
    end do
    do n=1,compact%nop(1)
      if( .not. compact%tfx(n)%null_op ) call compact%tfx(n)%remove()
    end do
    do n=1,compact%nop(3)
      if( .not. compact%sfz(n)%null_op ) call compact%sfz(n)%remove()
    end do
    do n=1,compact%nop(2)
      if( .not. compact%sfy(n)%null_op ) call compact%sfy(n)%remove()
    end do
    do n=1,compact%nop(1)
      if( .not. compact%sfx(n)%null_op ) call compact%sfx(n)%remove()
    end do
    do n=1,compact%nop(3)
      if( .not. compact%gfz(n)%null_op ) call compact%gfz(n)%remove()
    end do
    do n=1,compact%nop(2)
      if( .not. compact%gfy(n)%null_op ) call compact%gfy(n)%remove()
    end do
    do n=1,compact%nop(1)
      if( .not. compact%gfx(n)%null_op ) call compact%gfx(n)%remove()
    end do
    do n=1,compact%nop(3)
      if( .not. compact%d8z(n)%null_op ) call compact%d8z(n)%remove()
    end do
    do n=1,compact%nop(2)
      if( .not. compact%d8y(n)%null_op ) call compact%d8y(n)%remove()
    end do
    do n=1,compact%nop(1)
      if( .not. compact%d8x(n)%null_op ) call compact%d8x(n)%remove()
    end do
    do n=1,compact%nop(3)
      if( .not. compact%d4z(n)%null_op ) call compact%d4z(n)%remove()
    end do
    do n=1,compact%nop(2)
      if( .not. compact%d4y(n)%null_op ) call compact%d4y(n)%remove()
    end do
    do n=1,compact%nop(1)
      if( .not. compact%d4x(n)%null_op ) call compact%d4x(n)%remove()
    end do
    do n=1,compact%nop(3)
      if( .not. compact%d2z(n)%null_op ) call compact%d2z(n)%remove()
    end do
    do n=1,compact%nop(2)
      if( .not. compact%d2y(n)%null_op ) call compact%d2y(n)%remove()
    end do
    do n=1,compact%nop(1)
      if( .not. compact%d2x(n)%null_op ) call compact%d2x(n)%remove()
    end do
    do n=1,compact%nop(3)
      if( .not. compact%d1z(n)%null_op ) call compact%d1z(n)%remove()
    end do
    do n=1,compact%nop(2)
      if( .not. compact%d1y(n)%null_op ) call compact%d1y(n)%remove()
    end do
    do n=1,compact%nop(1)
      if( .not. compact%d1x(n)%null_op ) call compact%d1x(n)%remove()
    end do
    compact%mbc = 0 ; compact%nop = 1
  end subroutine remove_compact_ops
  
   SUBROUTINE setup_ghost_buffers(patch) ! on base grid
    IMPLICIT NONE
     type(patch_type), intent(in) :: patch
     integer :: ax,ay,az,px,py,pz
     ax = patch%ax ; ay = patch%ay ; az = patch%az
     px = patch%px ; py = patch%py ; pz = patch%pz
     allocate(vbr1x(4,ay,az),vbr2x(4,ay,az),vbs1x(4,ay,az),vbs2x(4,ay,az),dvopx(4,ay,az),dvox(4,ay,az,0:px-1))
     allocate(vbr1y(ax,4,az),vbr2y(ax,4,az),vbs1y(ax,4,az),vbs2y(ax,4,az),dvopy(4,ax,az),dvoy(4,ax,az,0:py-1))
     allocate(vbr1z(ax,ay,4),vbr2z(ax,ay,4),vbs1z(ax,ay,4),vbs2z(ax,ay,4),dvopz(4,ax,ay),dvoz(4,ax,ay,0:pz-1))
     !$omp target enter data map(alloc:vbr1x,vbr2x,vbs1x,vbs2x,dvopx,dvox, &
     !$omp                             vbr1y,vbr2y,vbs1y,vbs2y,dvopy,dvoy, &
     !$omp                             vbr1z,vbr2z,vbs1z,vbs2z,dvopz,dvoz) if(gpu_kernel==1)
   END SUBROUTINE setup_ghost_buffers

   SUBROUTINE remove_ghost_buffers() ! on base grid
    IMPLICIT NONE     
     !$omp target exit data map(delete:vbr1x,vbr2x,vbs1x,vbs2x,dvopx,dvox, &
     !$omp                             vbr1y,vbr2y,vbs1y,vbs2y,dvopy,dvoy, &
     !$omp                             vbr1z,vbr2z,vbs1z,vbs2z,dvopz,dvoz) if(gpu_kernel==1)
     deallocate(vbr1z,vbr2z,vbs1z,vbs2z,dvopz,dvoz)
     deallocate(vbr1y,vbr2y,vbs1y,vbs2y,dvopy,dvoy)
     deallocate(vbr1x,vbr2x,vbs1x,vbs2x,dvopx,dvox)
   END SUBROUTINE remove_ghost_buffers

end module LES_compact

