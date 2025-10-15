module LES_compact_operators
  USE iso_c_binding
  !USE mapp_exosim_annotation, ONLY : exosim_annotation_begin,exosim_annotation_end
  !use les_input, only : gpu_kernel
  use LES_objects, only : compact_ops=>compact_ptr      ! , mesh_data=>mesh_ptr
										
  IMPLICIT NONE
  REAL(KIND=c_double), PARAMETER :: zero=0.0_c_double
  INTEGER :: gpu_kernel = 1
  
contains

  subroutine d1x(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    !call exosm_annotation_begin("d1x")
    !$omp target data map(to:v) map(from:dv) if(gpu_kernel==1)			! Okay to leave redundant data maps?
    ! perform operation?
    if( compact_ops%control%null_opx ) then
      !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
    	do i=1,ax; do j=1,ay; do k=1,az
      	dv(i,j,k) = zero
     	end do; end do; end do;
      !$omp end target teams distribute parallel do
    else
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. associated(compact_ops%d1x(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
   	call compact_ops%d1x(iop)%evalx(ax,ay,az,v,dv,vb1,vb2)

    ! apply metric here or later? here d is scalar or size(dv,1)
    !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
    do k=1,az; do j=1,ay
     dv(:,j,k) = dv(:,j,k)/compact_ops%dx
     ! dv = dv/mesh_data%d1  ! 3D metric
    enddo; enddo
		!$omp end target teams distribute parallel do
		end if
		!$omp end target data
		!call exosm_annotation_end("d1x")
  end subroutine d1x

  subroutine d1y(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    !call exosm_annotation_begin("d1y")
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    !$omp target data map(to:v) map(from:dv) if(gpu_kernel==1)			! Okay to leave redundant data maps?
    ! perform operation?
    if( compact_ops%control%null_opy ) then
      !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
    	do i=1,ax; do j=1,ay; do k=1,az
      	dv(i,j,k) = zero
     	end do; end do; end do;
      !$omp end target teams distribute parallel do
    else
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. associated(compact_ops%d1y(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    call compact_ops%d1y(iop)%evaly(ax,ay,az,v,dv,vb1,vb2)
    
    ! apply metric here or later? here d is scalar or size(dv,2)
    !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
    do k=1,az; do i=1,ax
     dv(i,:,k) = dv(i,:,k)/compact_ops%dy
     !  dv = dv/mesh_data%d2  ! 3D metric
    enddo; enddo
		!$omp end target teams distribute parallel do
    end if
		!$omp end target data
		!call exosm_annotation_end("d1y")
  end subroutine d1y

  subroutine d1z(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in)  :: v
    real(kind=c_double), dimension(:,:,:), intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    !call exosm_annotation_begin("d1z")
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    !$omp target data map(to:v) map(from:dv) if(gpu_kernel==1)			! Okay to leave redundant data maps?
    ! perform operation?
    if( compact_ops%control%null_opz ) then
      !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
    	do i=1,ax; do j=1,ay; do k=1,az
      	dv(i,j,k) = zero
     	end do; end do; end do;
      !$omp end target teams distribute parallel do
    else
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. associated(compact_ops%d1z(2)%ar) ) iop = 2
    endif
      ! calculate grid derivative sans metric
      call compact_ops%d1z(iop)%evalz(ax,ay,az,v,dv,vb1,vb2)
		  
		  ! apply metric here or later? here d is scalar or size(dv,3)
      !$omp target teams distribute parallel do collapse(2) if(gpu_kernel==1)
      do j=1,ay;	do i=1,ax
       dv(i,j,:) = dv(i,j,:)/compact_ops%dz
       ! dv = dv/mesh_data%d3  ! 3D metric
      enddo; enddo
		  !$omp end target teams distribute parallel do
    end if
		!$omp end target data
		!call exosm_annotation_end("d1z")
  end subroutine d1z

  subroutine d2x(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    ! perform operation?
    if( compact_ops%control%null_opx ) then
      dv = zero
      return
    endif
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. associated(compact_ops%d2x(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    call compact_ops%d2x(iop)%evalx(ax,ay,az,v,dv,vb1,vb2)
    ! apply metric here or later? here d is scalar or size(dv,1)
    forall(j=1:size(dv,2),k=1:size(dv,3)) dv(:,j,k) = dv(:,j,k)/compact_ops%dx**2
    ! dv = dv/mesh_data%d1**2  ! 3D metric
  end subroutine d2x

  subroutine d2y(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    ! perform operation?
    if( compact_ops%control%null_opy ) then
      dv = zero
      return
    endif
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. associated(compact_ops%d2y(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    call compact_ops%d2y(iop)%evaly(ax,ay,az,v,dv,vb1,vb2)
    ! apply metric here or later? here d is scalar or size(dv,2)
    forall(i=1:size(dv,1),k=1:size(dv,3)) dv(i,:,k) = dv(i,:,k)/compact_ops%dy**2
    ! dv = dv/mesh_data%d2**2  ! 3D metric
  end subroutine d2y

  subroutine d2z(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in)  :: v
    real(kind=c_double), dimension(:,:,:), intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    ! perform operation?
    if( compact_ops%control%null_opz ) then
      dv = zero
      return
    endif
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. associated(compact_ops%d2z(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    call compact_ops%d2z(iop)%evalz(ax,ay,az,v,dv,vb1,vb2)
    ! apply metric here or later? here d is scalar or size(dv,3)
    forall(i=1:size(dv,1),j=1:size(dv,2)) dv(i,j,:) = dv(i,j,:)/compact_ops%dz**2
    ! dv = dv/mesh_data%d3**2  ! 3D metric
  end subroutine d2z

  subroutine d4x(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    ! perform operation?
    if( compact_ops%control%null_opx ) then
      dv = zero
      return
    endif
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. associated(compact_ops%d4x(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    call compact_ops%d4x(iop)%evalx(ax,ay,az,v,dv,vb1,vb2)
    ! No metric... unity assumed
  end subroutine d4x

  subroutine d4y(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    ! perform operation?
    if( compact_ops%control%null_opy ) then
      dv = zero
      return
    endif
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. associated(compact_ops%d4y(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    call compact_ops%d4y(iop)%evaly(ax,ay,az,v,dv,vb1,vb2)
    ! No metric... unity assumed
  end subroutine d4y

  subroutine d4z(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in)  :: v
    real(kind=c_double), dimension(:,:,:), intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    ! perform operation?
    if( compact_ops%control%null_opz ) then
      dv = zero
      return
    endif
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. associated(compact_ops%d4z(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    call compact_ops%d4z(iop)%evalz(ax,ay,az,v,dv,vb1,vb2)
    ! No metric... unity assumed
  end subroutine d4z

  subroutine d8x(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    !call exosm_annotation_begin("d8x")
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    !$omp target data map(to:v) map(from:dv) if(gpu_kernel==1)
    ! perform operation?
    if( compact_ops%control%null_opx ) then
      !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
    	do i=1,ax; do j=1,ay; do k=1,az
      	dv(i,j,k) = zero
     	end do; end do; end do;
      !$omp end target teams distribute parallel do
    else
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. associated(compact_ops%d8x(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    call compact_ops%d8x(iop)%evalx(ax,ay,az,v,dv,vb1,vb2)
    
    ! apply metric here or later? here d is scalar or size(dv,1)
 !   forall(j=1:size(dv,2),k=1:size(dv,3)) dv(:,j,k) = dv(:,j,k)/compact_ops%dx**8
    ! dv = dv/mesh_data%d1**8  ! 3D metric
    end if
    !$omp end target data
    !call exosm_annotation_end("d8x")
  end subroutine d8x

  subroutine d8y(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    !call exosm_annotation_begin("d8y")
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    !$omp target data map(to:v) map(from:dv) if(gpu_kernel==1)
    ! perform operation?
    if( compact_ops%control%null_opy ) then
    	!$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
    	do i=1,ax; do j=1,ay; do k=1,az
      	dv(i,j,k) = zero
     	end do; end do; end do;
      !$omp end target teams distribute parallel do
    else
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. associated(compact_ops%d8y(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    call compact_ops%d8y(iop)%evaly(ax,ay,az,v,dv,vb1,vb2)
    
    ! apply metric here or later? here d is scalar or size(dv,2)
 !   forall(i=1:size(dv,1),k=1:size(dv,3)) dv(i,:,k) = dv(i,:,k)/compact_ops%dy**8
    ! dv = dv/mesh_data%d2**8  ! 3D metric
    end if
    !$omp end target data
    !call exosm_annotation_end("d8y")
  end subroutine d8y

  subroutine d8z(v,dv,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in)  :: v
    real(kind=c_double), dimension(:,:,:), intent(out) :: dv
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int), intent(in), optional :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    !call exosm_annotation_begin("d8z")
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    !$omp target data map(to:v) map(from:dv) if(gpu_kernel==1)
    ! perform operation?
    if( compact_ops%control%null_opz ) then
      !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
    	do i=1,ax; do j=1,ay; do k=1,az
      	dv(i,j,k) = zero
     	end do; end do; end do;
      !$omp end target teams distribute parallel do
    else
    ! symmetry option
    iop = 1
    if( present(bc) ) then
      if( bc > 0 ) iop = bc
      if( bc == -1 .and. associated(compact_ops%d8z(2)%ar) ) iop = 2
    endif
    ! calculate grid derivative sans metric
    call compact_ops%d8z(iop)%evalz(ax,ay,az,v,dv,vb1,vb2)
    
    ! apply metric here or later? here d is scalar or size(dv,3)
!    forall(i=1:size(dv,1),j=1:size(dv,2)) dv(i,j,:) = dv(i,j,:)/compact_ops%dz**8
    ! dv = dv/mesh_data%d3**8  ! 3D metric
    end if
    !$omp end target data
    !call exosm_annotation_end("d8z")
  end subroutine d8z

  subroutine filterx(v,fv,nf,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: fv
    integer(c_int), intent(in) :: nf
    integer(c_int), intent(in), optional :: bc
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int) :: iop,i,j,k,ax,ay,az
    !call exosm_annotation_begin("filterx")
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    !$omp target data map(to:v) map(from:fv) if(gpu_kernel==1)
    ! perform operation?
    if( compact_ops%control%null_opx ) then
      !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
    	do i=1,ax; do j=1,ay; do k=1,az
        fv(i,j,k) = v(i,j,k)
      end do; end do; end do;
      !$omp end target teams distribute parallel do
    else
      ! symmetry option
      iop = 1
      if( present(bc) ) then
        if( bc > 0 ) iop = bc
        if( bc == -1 .and. compact_ops%nop(1) > 1 ) iop = 2
      endif
      ! select filter
      if( nf == compact_ops%control%gfspec ) then
        call compact_ops%gfx(iop)%evalx(ax,ay,az,v,fv,vb1,vb2)
      else if( nf == compact_ops%control%sfspec ) then
        call compact_ops%sfx(iop)%evalx(ax,ay,az,v,fv,vb1,vb2)
      else if( nf == compact_ops%control%tfspec ) then
        call compact_ops%tfx(iop)%evalx(ax,ay,az,v,fv,vb1,vb2)
      end if
    end if
    !$omp end target data
    !call exosm_annotation_end("filterx")
  end subroutine filterx

  subroutine filtery(v,fv,nf,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: fv
    integer(c_int), intent(in) :: nf
    integer(c_int), intent(in), optional :: bc
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int) :: iop,i,j,k,ax,ay,az
    !call exosm_annotation_begin("filtery")
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    !$omp target data map(to:v) map(from:fv) if(gpu_kernel==1)
    ! perform operation?
    if( compact_ops%control%null_opy ) then
      !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
    	do i=1,ax; do j=1,ay; do k=1,az
        fv(i,j,k) = v(i,j,k)
      end do; end do; end do;
      !$omp end target teams distribute parallel do
    else
      ! symmetry option
      iop = 1
      if( present(bc) ) then
        if( bc > 0 ) iop = bc
        if( bc == -1 .and. compact_ops%nop(2) > 1 ) iop = 2
      endif
      ! select filter
      if( nf == compact_ops%control%gfspec ) then
        call compact_ops%gfy(iop)%evaly(ax,ay,az,v,fv,vb1,vb2)
      else if( nf == compact_ops%control%sfspec ) then
        call compact_ops%sfy(iop)%evaly(ax,ay,az,v,fv,vb1,vb2)
      else if( nf == compact_ops%control%tfspec ) then
        call compact_ops%tfy(iop)%evaly(ax,ay,az,v,fv,vb1,vb2)
      end if
    endif
    !$omp end target data
    !call exosm_annotation_end("filtery")
  end subroutine filtery

  subroutine filterz(v,fv,nf,bc,vb1,vb2)
    IMPLICIT NONE
    real(kind=c_double), dimension(:,:,:), intent(in) :: v
    real(kind=c_double), dimension(:,:,:),intent(out) :: fv
    integer(c_int), intent(in) :: nf
    integer(c_int), intent(in), optional :: bc
    real(kind=c_double), dimension(:,:,:), intent(in), optional :: vb1,vb2
    integer(c_int) :: iop,i,j,k,ax,ay,az
    !call exosm_annotation_begin("filterz")
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    !$omp target data map(to:v) map(from:fv) if(gpu_kernel==1)
    ! perform operation?
    if( compact_ops%control%null_opz ) then
      !$omp target teams distribute parallel do collapse(3) if(gpu_kernel==1)
    	do i=1,ax; do j=1,ay; do k=1,az
        fv(i,j,k) = v(i,j,k)
      end do; end do; end do;
      !$omp end target teams distribute parallel do
    else
      ! symmetry option
      iop = 1
      if( present(bc) ) then
        if( bc > 0 ) iop = bc
        if( bc == -1 .and. compact_ops%nop(3) > 1 ) iop = 2
      endif
      ! select filter
      if( nf == compact_ops%control%gfspec ) then
        call compact_ops%gfz(iop)%evalz(ax,ay,az,v,fv,vb1,vb2)
      else if( nf == compact_ops%control%sfspec ) then
        call compact_ops%sfz(iop)%evalz(ax,ay,az,v,fv,vb1,vb2)
      else if( nf == compact_ops%control%tfspec ) then
        call compact_ops%tfz(iop)%evalz(ax,ay,az,v,fv,vb1,vb2)
      end if
    endif
    !$omp end target data
    !call exosm_annotation_end("filterz")
  end subroutine filterz



  SUBROUTINE icfx(v, iv, bc, vb1, vb2, iv1, iv2)
    IMPLICIT NONE
    REAL(kind=C_DOUBLE), DIMENSION(:, :, :), INTENT(in) :: v
    REAL(kind=C_DOUBLE), DIMENSION(:, :, :), INTENT(out) :: iv
    REAL(kind=C_DOUBLE), DIMENSION(:, :, :), INTENT(in), OPTIONAL :: vb1, vb2 ! ghost values
    REAL(kind=C_DOUBLE), DIMENSION(:, :), INTENT(in), OPTIONAL :: iv1, iv2 ! boundary values
    INTEGER(C_INT), INTENT(in), OPTIONAL :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    ! perform operation?
    IF (compact_ops%control%null_opx) THEN
!      iv = zero  ! ?
       RETURN
    END IF
    ! symmetry option
    iop = 1
    IF (PRESENT(bc)) THEN
       IF (bc > 0) iop = bc
       IF (bc == -1 .AND. ASSOCIATED(compact_ops%icfx(2)%ar)) iop = 2
    END IF
    call compact_ops%icfx(iop)%evalx(ax,ay,az,v, iv, vb1, vb2, iv1, iv2)
 END SUBROUTINE icfx

 SUBROUTINE icfy(v, iv, bc, vb1, vb2, iv1, iv2)
    IMPLICIT NONE
    REAL(kind=C_DOUBLE), DIMENSION(:, :, :), INTENT(in) :: v
    REAL(kind=C_DOUBLE), DIMENSION(:, :, :), INTENT(out) :: iv
    REAL(kind=C_DOUBLE), DIMENSION(:, :, :), INTENT(in), OPTIONAL :: vb1, vb2 ! ghost values
    REAL(kind=C_DOUBLE), DIMENSION(:, :), INTENT(in), OPTIONAL :: iv1, iv2 ! boundary values
    INTEGER(C_INT), INTENT(in), OPTIONAL :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    ! perform operation?
    IF (compact_ops%control%null_opy) THEN
!      iv = zero  ! ?
       RETURN
    END IF
    ! symmetry option
    iop = 1
    IF (PRESENT(bc)) THEN
       IF (bc > 0) iop = bc
       IF (bc == -1 .AND. ASSOCIATED(compact_ops%icfy(2)%ar)) iop = 2
    END IF
    call compact_ops%icfy(iop)%evaly(ax,ay,az,v, iv, vb1, vb2, iv1, iv2)
 END SUBROUTINE icfy

 SUBROUTINE icfz(v, iv, bc, vb1, vb2, iv1, iv2)
    IMPLICIT NONE
    REAL(kind=C_DOUBLE), DIMENSION(:, :, :), INTENT(in) :: v
    REAL(kind=C_DOUBLE), DIMENSION(:, :, :), INTENT(out) :: iv
    REAL(kind=C_DOUBLE), DIMENSION(:, :, :), INTENT(in), OPTIONAL :: vb1, vb2 ! ghost values
    REAL(kind=C_DOUBLE), DIMENSION(:, :), INTENT(in), OPTIONAL :: iv1, iv2 ! boundary values
    INTEGER(C_INT), INTENT(in), OPTIONAL :: bc
    integer(c_int) :: iop,i,j,k,ax,ay,az
    ax=size(v,1); ay=size(v,2); az=size(v,3)
    ! perform operation?
    IF (compact_ops%control%null_opz) THEN
!      iv = zero  ! ?
       RETURN
    END IF
    ! symmetry option
    iop = 1
    IF (PRESENT(bc)) THEN
       IF (bc > 0) iop = bc
       IF (bc == -1 .AND. ASSOCIATED(compact_ops%icfz(2)%ar)) iop = 2
    END IF
    call compact_ops%icfz(iop)%evalz(ax,ay,az,v, iv, vb1, vb2, iv1, iv2)
 END SUBROUTINE icfz




end module LES_compact_operators

