module LES_compact_r3
  USE iso_c_binding
  USE MPI
  !USE LES_input, ONLY : bpp_lus_opt,use_ppent_opt,directcom,gpu_kernel
  USE LES_stencils
  USE LES_comm,  ONLY : mpierr,mpistatus,master
  use LES_pentadiagonal, ONLY :   ppentlus, &
    bpentLUS3x, ppentLUS3x, bpentLUS3y, ppentLUS3y, bpentLUS3z, ppentLUS3z, &
    btrid_block4_lus, ptrid_block4_lus, bpentLUS2y 
  use LES_compact_basetype, only : compact_op1, &
  																 vbr1x,vbr2x,vbs1x,vbs2x,vbr1y,vbr2y,vbs1y,vbs2y, &
  																 vbr1z,vbr2z,vbs1z,vbs2z,dvopx,dvox,dvopy,dvoy, &
  																 dvopz,dvoz

  IMPLICIT NONE

  integer :: directcom = 1
  integer :: gpu_kernel = 1
  logical(c_bool) :: bpp_lus_opt = .true. , use_ppent_opt = .true.
  
  REAL(KIND=c_double), PARAMETER :: zero=0.0_c_double, one=1.0_c_double
  LOGICAL(c_bool) :: debug=.false.

  ! these extensions override the generic operations

  type, extends(compact_op1) :: compact_op1_r3  ! custom operator (backward compat.)
  contains
    procedure :: evalx => eval_compact_op1x_r3  ! rhs stencil = 3
    procedure :: evaly => eval_compact_op1y_r3  ! rhs stencil = 3
    procedure :: evalz => eval_compact_op1z_r3  ! rhs stencil = 3
  end type
  
contains

! "optimized" r3 operators for backward compatibility with matrix.f

  subroutine eval_compact_op1x_r3(op,ax,ay,az,v,dv,vb1,vb2,dv1,dv2)! nor=3, nol=2, uses 1D op type
    implicit none
    class(compact_op1_r3), intent(in) :: op
    integer, intent(in) :: ax,ay,az
    real(kind=c_double), dimension(ax,ay,az), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(ax,ay,az), intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    real(kind=c_double), dimension(op%ncr) :: vp
    real(kind=c_double), dimension(ax) :: sumr
    integer :: nb,nsr,i,j,k
    integer :: nor,nir,nr,nol,nl,ni,np  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v  ! filter / interp
!      print *,'null op in x'
      RETURN
    END IF
    !ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)
    !if( ax /= op%m ) then
    !  print *,'*** error: mismatch in x operation size ***',ax,op%m
    !  stop
    !endif
    if( op%nor /= 3 ) then
      print *,'*** error: mismatch in x stencil size ***',3,op%nor
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np
    ! explicit part
    allocate( vbr1(3,ay,az),vbr2(3,ay,az) )
! ghost data
   allocate( vbs1(3,ay,az),vbs2(3,ay,az) )
    if( np > 1 ) then  ! use parallel solver
      vbs2 = v(ax-2:ax,:,:)
      vbs1 = v(1:3,:,:)
      nsr = size(vbs1)
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
    else if( op%periodic ) then
      vbr1 = v(ax-2:ax,:,:)
      vbr2 = v(1:3,:,:)
      if( debug ) print *,'periodic in x'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,1)  ! assumes nb >= nor
        vbr1 = vb1(nb-2:nb,:,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,1)  ! assumes nb >= nor
        vbr2 = vb2(1:3,:,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    sumr = sum(op%ar,dim=1)
    if( op%lo == MPI_PROC_NULL .and. op%bc(1) == -1 ) then
    do k=1,az
    do j=1,ay
      dv(1,j,k) = sum(op%ar(1:3,1)*(vbr1(1:3,j,k)-v(1,j,k)))+sum(op%ar(4:7,1)*(v(1:4,j,k)-v(1,j,k)))+sumr(1)*v(1,j,k)
      dv(2,j,k) = sum(op%ar(1:2,2)*(vbr1(2:3,j,k)-v(2,j,k)))+sum(op%ar(3:7,2)*(v(1:5,j,k)-v(2,j,k)))+sumr(2)*v(2,j,k)
      dv(3,j,k) = sum(op%ar(1:1,3)*(vbr1(3:3,j,k)-v(3,j,k)))+sum(op%ar(2:7,3)*(v(1:6,j,k)-v(3,j,k)))+sumr(3)*v(3,j,k)
    end do
    end do
    else
    do k=1,az
    do j=1,ay
      dv(1,j,k) = sum(op%ar(1:3,1)*(vbr1(1:3,j,k)-v(1,j,k)))+sum(op%ar(4:7,1)*(v(1:4,j,k)-v(1,j,k)))
      dv(2,j,k) = sum(op%ar(1:2,2)*(vbr1(2:3,j,k)-v(2,j,k)))+sum(op%ar(3:7,2)*(v(1:5,j,k)-v(2,j,k)))
      dv(3,j,k) = sum(op%ar(1:1,3)*(vbr1(3:3,j,k)-v(3,j,k)))+sum(op%ar(2:7,3)*(v(1:6,j,k)-v(3,j,k)))
    end do
    end do
    endif
    do k=1,az
    do j=1,ay
      do i=4,ax-3
        vp = v(i-3:i+3,j,k)-v(i,j,k)
        dv(i,j,k) = sum(op%ar(:,i)*vp) ! +sumr(i)*v(i,j,k)
!        dv(i,j,k) = sum(op%ar(:,i)*v(i-3:i+3,j,k))
      end do
    end do
    end do
    if( op%hi == MPI_PROC_NULL .and. op%bc(2) == -1 ) then
    do k=1,az
    do j=1,ay
      dv(ax-2,j,k) = sum(op%ar(1:6,ax-2)*(v(ax-5:ax,j,k)-v(ax-2,j,k)))+sum(op%ar(7:7,ax-2)*(vbr2(1:1,j,k)-v(ax-2,j,k)))+sumr(ax-2)*v(ax-2,j,k)
      dv(ax-1,j,k) = sum(op%ar(1:5,ax-1)*(v(ax-4:ax,j,k)-v(ax-1,j,k)))+sum(op%ar(6:7,ax-1)*(vbr2(1:2,j,k)-v(ax-1,j,k)))+sumr(ax-1)*v(ax-1,j,k)
      dv(ax  ,j,k) = sum(op%ar(1:4,ax  )*(v(ax-3:ax,j,k)-v(ax  ,j,k)))+sum(op%ar(5:7,ax  )*(vbr2(1:3,j,k)-v(ax  ,j,k)))+sumr(ax  )*v(ax  ,j,k)
    end do
    end do
    else
    do k=1,az
    do j=1,ay
      dv(ax-2,j,k) = sum(op%ar(1:6,ax-2)*(v(ax-5:ax,j,k)-v(ax-2,j,k)))+sum(op%ar(7:7,ax-2)*(vbr2(1:1,j,k)-v(ax-2,j,k)))
      dv(ax-1,j,k) = sum(op%ar(1:5,ax-1)*(v(ax-4:ax,j,k)-v(ax-1,j,k)))+sum(op%ar(6:7,ax-1)*(vbr2(1:2,j,k)-v(ax-1,j,k)))
      dv(ax  ,j,k) = sum(op%ar(1:4,ax  )*(v(ax-3:ax,j,k)-v(ax  ,j,k)))+sum(op%ar(5:7,ax  )*(vbr2(1:3,j,k)-v(ax  ,j,k)))
    end do
    end do
    endif
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(1,:,:)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(ax,:,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) then
      if( op%null_option == 1 ) dv=dv+v  ! filter / interp
      return
    endif
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3x(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
    else
      call bpentLUS3x(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( np == 1 ) then
      if( op%null_option == 1 ) dv=dv+v  ! filter / interp
      return
    endif
    ! parallel solver
      allocate( dvop(4,ay,az), dvo(4,ay,az,0:np-1) )
      forall(i=1:2,j=1:ay,k=1:az) 
        dvop(i,j,k) = dv(i,j,k)
        dvop(i+2,j,k) = dv(ax+i-2,j,k)
      end forall
      if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
      if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
      nsr = size(dvop)
      select case( op%directcom )
      case( 1 ) ! mpi_allgather
        call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvo, np, ay, az )
        else
          call btrid_block4_lus( op%aa, dvo, np, ay, az )
        endif
        if( op%lo /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,1:2)*dvo(3:4,j,k,op%lo))
        if( op%hi /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,3:4)*dvo(1:2,j,k,op%hi))
      case( 2 ) ! mpi_gather/scatter
        call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvo, np, ay, az )
          else
            call btrid_block4_lus( op%aa, dvo, np, ay, az )
          endif
        else
          dvo = zero
        endif
        ! shuffle solution
        dvop(3:4,:,:) = dvo(1:2,:,:,0)
        do i=0,np-2
          dvop(1:2,:,:) = dvo(3:4,:,:,i)
          dvo(3:4,:,:,i) = dvo(1:2,:,:,i+1)
          dvo(1:2,:,:,i+1) = dvop(1:2,:,:)
        end do
        dvo(1:2,:,:,0) = dvo(3:4,:,:,np-1)
        dvo(3:4,:,:,np-1) = dvop(3:4,:,:)
        call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(i,:)*dvop(:,j,k))
      end select
      deallocate( dvop, dvo )
      if( op%null_option == 1 ) dv=dv+v  ! filter / interp
  end subroutine eval_compact_op1x_r3

  subroutine eval_compact_op1y_r3(op,ax,ay,az,v,dv,vb1,vb2,dv1,dv2) ! nor=3, nol=2, uses 1D op type
    implicit none
    class(compact_op1_r3), intent(in) :: op
    integer, intent(in) :: ax,ay,az
    real(kind=c_double), dimension(ax,ay,az), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(ax,ay,az), intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    real(kind=c_double), dimension(op%ncr) :: vp
    real(kind=c_double), dimension(ay) :: sumr
    integer :: nb,nsr,i,j,k
    integer :: nor,nir,nr,nol,nl,ni,np  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v  ! filter / interp
!      print *,'null op in y'
      RETURN
    END IF
    !ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)
    if( ay /= op%m ) then
      print *,'*** error: mismatch in y operation size ***',ay,op%m
      stop
    endif
    if( op%nor /= 3 ) then
      print *,'*** error: mismatch in y stencil size ***',3,op%nor
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np
    ! explicit part
    allocate( vbr1(ax,3,az),vbr2(ax,3,az) )
! ghost data
   allocate( vbs1(ax,3,az),vbs2(ax,3,az) )
    if( np > 1 ) then  ! use parallel solver
      vbs2 = v(:,ay-2:ay,:)
      vbs1 = v(:,1:3,:)
      nsr = size(vbs1)
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
    else if( op%periodic ) then
      vbr1 = v(:,ay-2:ay,:)
      vbr2 = v(:,1:3,:)
      if( debug ) print *,'periodic in y'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,2)  ! assumes nb >= nor
        vbr1 = vb1(:,nb-2:nb,:)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,2)  ! assumes nb >= nor
        vbr2 = vb2(:,1:3,:)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    sumr = sum(op%ar,dim=1)
    if( op%lo == MPI_PROC_NULL .and. op%bc(1) == -1 ) then
    do k=1,az
    do i=1,ax
      dv(i,1,k) = sum(op%ar(1:3,1)*(vbr1(i,1:3,k)-v(i,1,k)))+sum(op%ar(4:7,1)*(v(i,1:4,k)-v(i,1,k)))+sumr(1)*v(i,1,k)
      dv(i,2,k) = sum(op%ar(1:2,2)*(vbr1(i,2:3,k)-v(i,2,k)))+sum(op%ar(3:7,2)*(v(i,1:5,k)-v(i,2,k)))+sumr(2)*v(i,2,k)
      dv(i,3,k) = sum(op%ar(1:1,3)*(vbr1(i,3:3,k)-v(i,3,k)))+sum(op%ar(2:7,3)*(v(i,1:6,k)-v(i,3,k)))+sumr(3)*v(i,3,k)
    end do
    end do
    else
    do k=1,az
    do i=1,ax
      dv(i,1,k) = sum(op%ar(1:3,1)*(vbr1(i,1:3,k)-v(i,1,k)))+sum(op%ar(4:7,1)*(v(i,1:4,k)-v(i,1,k)))
      dv(i,2,k) = sum(op%ar(1:2,2)*(vbr1(i,2:3,k)-v(i,2,k)))+sum(op%ar(3:7,2)*(v(i,1:5,k)-v(i,2,k)))
      dv(i,3,k) = sum(op%ar(1:1,3)*(vbr1(i,3:3,k)-v(i,3,k)))+sum(op%ar(2:7,3)*(v(i,1:6,k)-v(i,3,k)))
    end do
    end do
    endif
    do k=1,az
    do j=4,ay-3
    do i=1,ax
        vp = v(i,j-3:j+3,k)-v(i,j,k)
        dv(i,j,k) = sum(op%ar(:,j)*vp) ! +sumr(j)*v(i,j,k)
!        dv(i,j,k) = sum(op%ar(:,j)*v(i,j-3:j+3,k))
    end do
    end do
    end do
    if( op%hi == MPI_PROC_NULL .and. op%bc(2) == -1 ) then
    do k=1,az
    do i=1,ax
      dv(i,ay-2,k) = sum(op%ar(1:6,ay-2)*(v(i,ay-5:ay,k)-v(i,ay-2,k)))+sum(op%ar(7:7,ay-2)*(vbr2(i,1:1,k)-v(i,ay-2,k)))+sumr(ay-2)*v(i,ay-2,k)
      dv(i,ay-1,k) = sum(op%ar(1:5,ay-1)*(v(i,ay-4:ay,k)-v(i,ay-1,k)))+sum(op%ar(6:7,ay-1)*(vbr2(i,1:2,k)-v(i,ay-1,k)))+sumr(ay-1)*v(i,ay-1,k)
      dv(i,ay  ,k) = sum(op%ar(1:4,ay  )*(v(i,ay-3:ay,k)-v(i,ay  ,k)))+sum(op%ar(5:7,ay  )*(vbr2(i,1:3,k)-v(i,ay  ,k)))+sumr(ay  )*v(i,ay  ,k)
    end do
    end do
    else
    do k=1,az
    do i=1,ax
      dv(i,ay-2,k) = sum(op%ar(1:6,ay-2)*(v(i,ay-5:ay,k)-v(i,ay-2,k)))+sum(op%ar(7:7,ay-2)*(vbr2(i,1:1,k)-v(i,ay-2,k)))
      dv(i,ay-1,k) = sum(op%ar(1:5,ay-1)*(v(i,ay-4:ay,k)-v(i,ay-1,k)))+sum(op%ar(6:7,ay-1)*(vbr2(i,1:2,k)-v(i,ay-1,k)))
      dv(i,ay  ,k) = sum(op%ar(1:4,ay  )*(v(i,ay-3:ay,k)-v(i,ay  ,k)))+sum(op%ar(5:7,ay  )*(vbr2(i,1:3,k)-v(i,ay  ,k)))
    end do
    end do
    endif
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,1,:)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,ay,:)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) then
      if( op%null_option == 1 ) dv=dv+v  ! filter / interp
      return
    endif
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3y(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
    else
      call bpentLUS3y(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( np == 1 ) then
      if( op%null_option == 1 ) dv=dv+v  ! filter / interp
      return
    endif
    ! parallel solver
      allocate( dvop(4,ax,az), dvo(4,ax,az,0:np-1) )
      forall(i=1:ax,j=1:2,k=1:az) 
        dvop(j,i,k) = dv(i,j,k)
        dvop(j+2,i,k) = dv(i,ay+j-2,k)
      end forall
      if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
      if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
      nsr = size(dvop)
      select case( op%directcom )
      case( 1 ) ! mpi_allgather
        call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvo, np, ax, az )
        else
          call btrid_block4_lus( op%aa, dvo, np, ax, az )
        endif
        if( op%lo /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,1:2)*dvo(3:4,i,k,op%lo))
        if( op%hi /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,3:4)*dvo(1:2,i,k,op%hi))
      case( 2 ) ! mpi_gather/scatter
        call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvo, np, ax, az )
          else
            call btrid_block4_lus( op%aa, dvo, np, ax, az )
          endif
        else
          dvo = zero
        endif
        ! shuffle solution
        dvop(3:4,:,:) = dvo(1:2,:,:,0)
        do i=0,np-2
          dvop(1:2,:,:) = dvo(3:4,:,:,i)
          dvo(3:4,:,:,i) = dvo(1:2,:,:,i+1)
          dvo(1:2,:,:,i+1) = dvop(1:2,:,:)
        end do
        dvo(1:2,:,:,0) = dvo(3:4,:,:,np-1)
        dvo(3:4,:,:,np-1) = dvop(3:4,:,:)
        call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(j,:)*dvop(:,i,k))
      end select
      deallocate( dvop, dvo )
      if( op%null_option == 1 ) dv=dv+v  ! filter / interp
  end subroutine eval_compact_op1y_r3

  subroutine eval_compact_op1z_r3(op,ax,ay,az,v,dv,vb1,vb2,dv1,dv2) ! nor=3, nol=2, uses 1D op type
    implicit none
    class(compact_op1_r3), intent(in) :: op
    integer, intent(in) :: ax,ay,az
    real(kind=c_double), dimension(ax,ay,az), intent(in) :: v
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2 ! ghost values
    real(kind=c_double), dimension(:,:), intent(in), optional :: dv1,dv2 ! boundary values
    real(kind=c_double), dimension(ax,ay,az), intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), allocatable :: vbr1,vbr2,vbs1,vbs2
    real(kind=c_double), dimension(:,:,:), allocatable :: dvop
    real(kind=c_double), dimension(:,:,:,:), allocatable :: dvo
    real(kind=c_double), dimension(op%ncr) :: vp
    real(kind=c_double), dimension(az) :: sumr
    integer :: nb,nsr,i,j,k
    integer :: nor,nir,nr,nol,nl,ni,np  ! surrogates
    IF (op%null_op) THEN
      if( op%null_option == 0 ) dv = zero
      if( op%null_option == 1 ) dv = v  ! filter / interp
!      print *,'null op in z'
      RETURN
    END IF
    !ax = size(v,1) ; ay = size(v,2) ; az = size(v,3)
    if( az /= op%m ) then
      print *,'*** error: mismatch in z operation size ***',az,op%m
      stop
    endif
    if( op%nor /= 3 ) then
      print *,'*** error: mismatch in z stencil size ***',3,op%nor
      stop
    endif
    nor = op%nor ; nir = op%nir ; nr = op%ncr
    nol = op%nol ; nl = op%ncl ; ni = op%nci
    np =  op%np
    ! explicit part
! ghost data
    allocate( vbr1(ax,ay,3),vbr2(ax,ay,3) )
    allocate( vbs1(ax,ay,3),vbs2(ax,ay,3) )
    if( np > 1 ) then  ! use parallel solver
      vbs2 = v(:,:,az-2:az)
      vbs1 = v(:,:,1:3)
      nsr = size(vbs1)
      call MPI_Sendrecv( vbs2, nsr, MPI_DOUBLE_PRECISION, op%hi, 0, &
                         vbr1, nsr, MPI_DOUBLE_PRECISION, op%lo, 0, &
                         op%hash, mpistatus, mpierr )
      call MPI_Sendrecv( vbs1, nsr, MPI_DOUBLE_PRECISION, op%lo, 1, &
                         vbr2, nsr, MPI_DOUBLE_PRECISION, op%hi, 1, &
                         op%hash, mpistatus, mpierr )
    else if( op%periodic ) then
      vbr1 = v(:,:,az-2:az)
      vbr2 = v(:,:,1:3)
      if( debug ) print *,'periodic in z'
    endif
    if( op%lo == MPI_PROC_NULL ) then ! lower non-periodic physical boundary
      if( present(vb1) ) then  ! use externally supplied values  ! .and. op%bc(1) = 2
        nb = size(vb1,3)  ! assumes nb >= nor
        vbr1 = vb1(:,:,nb-2:nb)
      else
        vbr1 = zero  ! no data or symmetry
      endif
    endif
    if( op%hi == MPI_PROC_NULL ) then  ! upper non-periodic physical boundary
      if( present(vb2) ) then  ! use externally supplied values  ! .and. op%bc(2) = 2
        nb = size(vb2,3)  ! assumes nb >= nor
        vbr2 = vb2(:,:,1:3)
      else
        vbr2 = zero  ! no data or symmetry
      endif
    endif
    deallocate( vbs1,vbs2 )
    sumr = sum(op%ar,dim=1)
    if( op%lo == MPI_PROC_NULL .and. op%bc(1) == -1 ) then
    do j=1,ay
    do i=1,ax
      dv(i,j,1) = sum(op%ar(1:3,1)*(vbr1(i,j,1:3)-v(i,j,1)))+sum(op%ar(4:7,1)*(v(i,j,1:4)-v(i,j,1)))+sumr(1)*v(i,j,1)
      dv(i,j,2) = sum(op%ar(1:2,2)*(vbr1(i,j,2:3)-v(i,j,2)))+sum(op%ar(3:7,2)*(v(i,j,1:5)-v(i,j,2)))+sumr(2)*v(i,j,2)
      dv(i,j,3) = sum(op%ar(1:1,3)*(vbr1(i,j,3:3)-v(i,j,3)))+sum(op%ar(2:7,3)*(v(i,j,1:6)-v(i,j,3)))+sumr(3)*v(i,j,3)
    end do
    end do
    else
    do j=1,ay
    do i=1,ax
      dv(i,j,1) = sum(op%ar(1:3,1)*(vbr1(i,j,1:3)-v(i,j,1)))+sum(op%ar(4:7,1)*(v(i,j,1:4)-v(i,j,1)))
      dv(i,j,2) = sum(op%ar(1:2,2)*(vbr1(i,j,2:3)-v(i,j,2)))+sum(op%ar(3:7,2)*(v(i,j,1:5)-v(i,j,2)))
      dv(i,j,3) = sum(op%ar(1:1,3)*(vbr1(i,j,3:3)-v(i,j,3)))+sum(op%ar(2:7,3)*(v(i,j,1:6)-v(i,j,3)))
    end do
    end do
    endif
    do k=4,az-3
    do j=1,ay
    do i=1,ax
        vp = v(i,j,k-3:k+3)-v(i,j,k)
        dv(i,j,k) = sum(op%ar(:,k)*vp) ! +sumr(k)*v(i,j,k)
!        dv(i,j,k) = sum(op%ar(:,k)*v(i,j,k-3:k+3))
    end do
    end do
    end do
    if( op%hi == MPI_PROC_NULL .and. op%bc(2) == -1 ) then
    do j=1,ay
    do i=1,ax
      dv(i,j,az-2) = sum(op%ar(1:6,az-2)*(v(i,j,az-5:az)-v(i,j,az-2)))+sum(op%ar(7:7,az-2)*(vbr2(i,j,1:1)-v(i,j,az-2)))+sumr(az-2)*v(i,j,az-2)
      dv(i,j,az-1) = sum(op%ar(1:5,az-1)*(v(i,j,az-4:az)-v(i,j,az-1)))+sum(op%ar(6:7,az-1)*(vbr2(i,j,1:2)-v(i,j,az-1)))+sumr(az-1)*v(i,j,az-1)
      dv(i,j,az  ) = sum(op%ar(1:4,az  )*(v(i,j,az-3:az)-v(i,j,az  )))+sum(op%ar(5:7,az  )*(vbr2(i,j,1:3)-v(i,j,az  )))+sumr(az  )*v(i,j,az  )
    end do
    end do
    else
    do j=1,ay
    do i=1,ax
      dv(i,j,az-2) = sum(op%ar(1:6,az-2)*(v(i,j,az-5:az)-v(i,j,az-2)))+sum(op%ar(7:7,az-2)*(vbr2(i,j,1:1)-v(i,j,az-2)))
      dv(i,j,az-1) = sum(op%ar(1:5,az-1)*(v(i,j,az-4:az)-v(i,j,az-1)))+sum(op%ar(6:7,az-1)*(vbr2(i,j,1:2)-v(i,j,az-1)))
      dv(i,j,az  ) = sum(op%ar(1:4,az  )*(v(i,j,az-3:az)-v(i,j,az  )))+sum(op%ar(5:7,az  )*(vbr2(i,j,1:3)-v(i,j,az  )))
    end do
    end do
    endif
    deallocate( vbr1,vbr2 )
    ! this is only appropriate for explicit bc
    if( (op%lo == MPI_PROC_NULL) .and. abs(op%bc(1)) /= 1 .and. present(dv1)) dv(:,:,1)=dv1   ! supply lower solution
    if( (op%hi == MPI_PROC_NULL) .and. abs(op%bc(2)) /= 1 .and. present(dv2)) dv(:,:,az)=dv2  ! supply upper solution
    if( .not. op%implicit_op ) then
      if( op%null_option == 1 ) dv=dv+v  ! filter / interp
      return
    endif
    ! implicit part
    if( np == 1 .and. op%periodic ) then
      call ppentLUS3z(op%al,dv,op%m,ax,ay,az)  ! periodic solution on a single process
    else
      call bpentLUS3z(op%al,dv,op%m,ax,ay,az)  ! locally non-periodic solution
    endif
    if( np == 1 ) then
      if( op%null_option == 1 ) dv=dv+v  ! filter / interp
      return
    endif
    ! parallel solver
      allocate( dvop(4,ax,ay), dvo(4,ax,ay,0:np-1) )
      forall(i=1:ax,j=1:ay,k=1:2) 
        dvop(k,i,j) = dv(i,j,k)
        dvop(k+2,i,j) = dv(i,j,az+k-2)
      end forall
      if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
      if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
      nsr = size(dvop)
      select case( op%directcom )
      case( 1 ) ! mpi_allgather
        call mpi_allgather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,op%hash,mpierr)
        if( op%periodic ) then
          call ptrid_block4_lus( op%aa, dvo, np, ax, ay )
        else
          call btrid_block4_lus( op%aa, dvo, np, ax, ay )
        endif
        if( op%lo /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,1:2)*dvo(3:4,i,j,op%lo))
        if( op%hi /= MPI_PROC_NULL ) forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,3:4)*dvo(1:2,i,j,op%hi))
      case( 2 ) ! mpi_gather/scatter
        call mpi_gather(dvop,nsr,MPI_DOUBLE_PRECISION,dvo,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%id == 0 ) then  ! only master solves
          if( op%periodic ) then
            call ptrid_block4_lus( op%aa, dvo, np, ax, ay )
          else
            call btrid_block4_lus( op%aa, dvo, np, ax, ay )
          endif
        else
          dvo = zero
        endif
        ! shuffle solution
        dvop(3:4,:,:) = dvo(1:2,:,:,0)
        do i=0,np-2
          dvop(1:2,:,:) = dvo(3:4,:,:,i)
          dvo(3:4,:,:,i) = dvo(1:2,:,:,i+1)
          dvo(1:2,:,:,i+1) = dvop(1:2,:,:)
        end do
        dvo(1:2,:,:,0) = dvo(3:4,:,:,np-1)
        dvo(3:4,:,:,np-1) = dvop(3:4,:,:)
        call mpi_scatter(dvo,nsr,MPI_DOUBLE_PRECISION,dvop,nsr,MPI_DOUBLE_PRECISION,0,op%hash,mpierr)
        if( op%lo == MPI_PROC_NULL ) dvop(1:2,:,:) = zero
        if( op%hi == MPI_PROC_NULL ) dvop(3:4,:,:) = zero
        forall(i=1:ax,j=1:ay,k=1:az) dv(i,j,k) = dv(i,j,k)-sum(op%rc(k,:)*dvop(:,i,j))
      end select
      deallocate( dvop, dvo )
      if( op%null_option == 1 ) dv=dv+v  ! filter / interp
  end subroutine eval_compact_op1z_r3

end module LES_compact_r3
