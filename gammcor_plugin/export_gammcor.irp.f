program export_trexio
  implicit none
  read_wf = .True.
  SOFT_TOUCH read_wf
  call run
end

subroutine run
  use trexio
  implicit none
  BEGIN_DOC
  !     Exports the wave function in TREXIO format
  END_DOC

  integer(trexio_t)              :: f ! TREXIO file handle
  integer(trexio_exit_code)      :: rc
  double precision, allocatable  :: factor(:)

  print *, 'TREXIO file : '//trim(trexio_file)
  print *, ''

  call system('cp -f '//trim(trexio_file)//' '//trim(trexio_file)//'.bak')
  if (backend == 0) then
    f = trexio_open(trexio_file, 'u', TREXIO_HDF5, rc)
  else if (backend == 1) then
    f = trexio_open(trexio_file, 'u', TREXIO_TEXT, rc)
  endif
  if (f == 0_8) then
    print *, 'Unable to open TREXIO file for writing'
    print *, 'rc = ', rc
    stop -1
  endif

! ------------------------------------------------------------------------------

! Electrons
! ---------

  print *, 'Electrons'

  rc = trexio_write_electron_up_num(f, elec_alpha_num)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_electron_dn_num(f, elec_beta_num)
  call trexio_assert(rc, TREXIO_SUCCESS)


! Nuclei
! ------

  print *, 'Nuclei'

  rc = trexio_write_nucleus_num(f, nucl_num)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_nucleus_charge(f, nucl_charge)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_nucleus_coord(f, nucl_coord_transp)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_nucleus_label(f, nucl_label, 32)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_nucleus_repulsion(f, nuclear_repulsion)
  call trexio_assert(rc, TREXIO_SUCCESS)


! Pseudo-potentials
! -----------------

  if (do_pseudo) then

    print *, 'ECP'
    integer                       :: num

    num = 0
    do k=1,pseudo_klocmax
      do i=1,nucl_num
        if (pseudo_dz_k(i,k) /= 0.d0) then
          num = num+1
        end if
      end do
    end do

    do l=0,pseudo_lmax
      do k=1,pseudo_kmax
        do i=1,nucl_num
          if (pseudo_dz_kl(i,k,l) /= 0.d0) then
            num = num+1
          end if
        end do
      end do
    end do

    integer, allocatable          :: ang_mom(:), nucleus_index(:), power(:), lmax(:)
    double precision, allocatable :: exponent(:), coefficient(:)

    allocate(ang_mom(num), nucleus_index(num), exponent(num), coefficient(num), power(num), &
            lmax(nucl_num) )

    do i=1,nucl_num
      lmax(i) = -1
      do l=0,pseudo_lmax
        do k=1,pseudo_kmax
          if (pseudo_dz_kl_transp(k,l,i) /= 0.d0) then
            lmax(i) = max(lmax(i), l)
          end if
        end do
      end do
    end do

    j = 0
    do i=1,nucl_num
      do k=1,pseudo_klocmax
        if (pseudo_dz_k_transp(k,i) /= 0.d0) then
          j = j+1
          ang_mom(j) = lmax(i)+1
          nucleus_index(j) = i
          exponent(j) = pseudo_dz_k_transp(k,i)
          coefficient(j) = pseudo_v_k_transp(k,i)
          power(j) = pseudo_n_k_transp(k,i)
        end if
      end do

      do l=0,lmax(i)
        do k=1,pseudo_kmax
          if (pseudo_dz_kl_transp(k,l,i) /= 0.d0) then
            j = j+1
            ang_mom(j) = l
            nucleus_index(j) = i
            exponent(j) = pseudo_dz_kl_transp(k,l,i)
            coefficient(j) = pseudo_v_kl_transp(k,l,i)
            power(j) = pseudo_n_kl_transp(k,l,i)
          end if
        end do
      end do
    end do


    lmax(:) = lmax(:)+1
    rc = trexio_write_ecp_max_ang_mom_plus_1(f, lmax)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_z_core(f, int(nucl_charge_remove))
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_num(f, num)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_ang_mom(f, ang_mom)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_nucleus_index(f, nucleus_index)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_exponent(f, exponent)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_coefficient(f, coefficient)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ecp_power(f, power)
    call trexio_assert(rc, TREXIO_SUCCESS)

  endif


! Basis
! -----

  print *, 'Basis'


  rc = trexio_write_basis_type(f, 'Gaussian', len('Gaussian'))
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_basis_prim_num(f, prim_num)
  call trexio_assert(rc, TREXIO_SUCCESS)

   rc = trexio_write_basis_shell_num(f, shell_num)
   call trexio_assert(rc, TREXIO_SUCCESS)

   rc = trexio_write_basis_nucleus_index(f, basis_nucleus_index)
   call trexio_assert(rc, TREXIO_SUCCESS)

   rc = trexio_write_basis_shell_ang_mom(f, shell_ang_mom)
   call trexio_assert(rc, TREXIO_SUCCESS)

   allocate(factor(shell_num))
   if (ao_normalized) then
     factor(1:shell_num) = shell_normalization_factor(1:shell_num)
   else
     factor(1:shell_num) = 1.d0
   endif
   rc = trexio_write_basis_shell_factor(f, factor)
   call trexio_assert(rc, TREXIO_SUCCESS)

   deallocate(factor)

  rc = trexio_write_basis_shell_index(f, shell_index)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_basis_exponent(f, prim_expo)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_basis_coefficient(f, prim_coef)
  call trexio_assert(rc, TREXIO_SUCCESS)

  allocate(factor(prim_num))
  if (primitives_normalized) then
    factor(1:prim_num) = prim_normalization_factor(1:prim_num)
  else
    factor(1:prim_num) = 1.d0
  endif
  rc = trexio_write_basis_prim_factor(f, factor)
  call trexio_assert(rc, TREXIO_SUCCESS)
  deallocate(factor)


! Atomic orbitals
! ---------------

  print *, 'AOs'

  rc = trexio_write_ao_num(f, ao_num)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_ao_cartesian(f, 1)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_ao_shell(f, ao_shell)
  call trexio_assert(rc, TREXIO_SUCCESS)

  integer :: i, pow0(3), powA(3), j, k, l, nz
  double precision :: normA, norm0, C_A(3), overlap_x, overlap_z, overlap_y, c
  nz=100

  C_A(1) = 0.d0
  C_A(2) = 0.d0
  C_A(3) = 0.d0

  allocate(factor(ao_num))
  if (ao_normalized) then
    do i=1,ao_num
      l = ao_first_of_shell(ao_shell(i))
      factor(i) = (ao_coef_normalized(i,1)+tiny(1.d0))/(ao_coef_normalized(l,1)+tiny(1.d0))
    enddo
  else
    factor(:) = 1.d0
  endif
  rc = trexio_write_ao_normalization(f, factor)
  call trexio_assert(rc, TREXIO_SUCCESS)
  deallocate(factor)

! One-e AO integrals
! ------------------

  if (export_ao_one_e_ints) then
    print *, 'AO one-e integrals'

    rc = trexio_write_ao_1e_int_overlap(f,ao_overlap)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ao_1e_int_kinetic(f,ao_kinetic_integrals)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_ao_1e_int_potential_n_e(f,ao_integrals_n_e)
    call trexio_assert(rc, TREXIO_SUCCESS)

    if (do_pseudo) then
      rc = trexio_write_ao_1e_int_ecp(f, ao_pseudo_integrals_local + ao_pseudo_integrals_non_local)
      call trexio_assert(rc, TREXIO_SUCCESS)
    endif

    rc = trexio_write_ao_1e_int_core_hamiltonian(f,ao_one_e_integrals)
    call trexio_assert(rc, TREXIO_SUCCESS)
  end if

! Two-e AO integrals
! ------------------

  if (export_ao_two_e_ints) then
    print *, 'AO two-e integrals'
    PROVIDE ao_two_e_integrals_in_map

    integer(8), parameter :: BUFSIZE=10000_8
    double precision :: eri_buffer(BUFSIZE), integral
    integer(4) :: eri_index(4,BUFSIZE)
    integer(8) :: icount, offset

    double precision, external :: get_ao_two_e_integral


    icount = 0_8
    offset = 0_8
    do l=1,ao_num
      do k=1,ao_num
        do j=l,ao_num
          do i=k,ao_num
            if (i==j .and. k<l) cycle
            if (i<j) cycle
            integral = get_ao_two_e_integral(i,j,k,l,ao_integrals_map)
            if (integral == 0.d0) cycle
            icount += 1_8
            eri_buffer(icount) = integral
            eri_index(1,icount) = i
            eri_index(2,icount) = j
            eri_index(3,icount) = k
            eri_index(4,icount) = l
            if (icount == BUFSIZE) then
              rc = trexio_write_ao_2e_int_eri(f, offset, icount, eri_index, eri_buffer)
              call trexio_assert(rc, TREXIO_SUCCESS)
              offset += icount
              icount = 0_8
            end if
          end do
        end do
      end do
    end do

    if (icount >= 0_8) then
      rc = trexio_write_ao_2e_int_eri(f, offset, icount, eri_index, eri_buffer)
      call trexio_assert(rc, TREXIO_SUCCESS)
    end if
  end if


! Molecular orbitals
! ------------------

  print *, 'MOs'

  rc = trexio_write_mo_type(f, mo_label, len(trim(mo_label)))
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_mo_num(f, mo_num)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_mo_coefficient(f, mo_coef)
  call trexio_assert(rc, TREXIO_SUCCESS)

!  if (trim(mo_label) == 'Canonical') then
!    rc = trexio_write_mo_energy(f, fock_matrix_diag_mo)
!    call trexio_assert(rc, TREXIO_SUCCESS)
!  endif


! One-e MO integrals
! ------------------

  if (export_mo_two_e_ints) then
    print *, 'MO one-e integrals'

    rc = trexio_write_mo_1e_int_kinetic(f,mo_kinetic_integrals)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_write_mo_1e_int_potential_n_e(f,mo_integrals_n_e)
    call trexio_assert(rc, TREXIO_SUCCESS)

    if (do_pseudo) then
      rc = trexio_write_mo_1e_int_ecp(f,mo_pseudo_integrals_local)
      call trexio_assert(rc, TREXIO_SUCCESS)
    endif

    rc = trexio_write_mo_1e_int_core_hamiltonian(f,mo_one_e_integrals)
    call trexio_assert(rc, TREXIO_SUCCESS)
  end if

! Two-e MO integrals
! ------------------

  if (export_mo_two_e_ints) then
    print *, 'MO two-e integrals'
    PROVIDE mo_two_e_integrals_in_map

    double precision, external :: mo_two_e_integral


    icount = 0_8
    offset = 0_8
    do l=1,mo_num
      do k=1,mo_num
        do j=l,mo_num
          do i=k,mo_num
            if (i==j .and. k<l) cycle
            if (i<j) cycle
            integral = mo_two_e_integral(i,j,k,l)
            if (integral == 0.d0) cycle
            icount += 1_8
            eri_buffer(icount) = integral
            eri_index(1,icount) = i
            eri_index(2,icount) = j
            eri_index(3,icount) = k
            eri_index(4,icount) = l
            if (icount == BUFSIZE) then
              rc = trexio_write_mo_2e_int_eri(f, offset, icount, eri_index, eri_buffer)
              call trexio_assert(rc, TREXIO_SUCCESS)
              offset += icount
              icount = 0_8
            end if
          end do
        end do
      end do
    end do

    if (icount > 0_8) then
      rc = trexio_write_mo_2e_int_eri(f, offset, icount, eri_index, eri_buffer)
      call trexio_assert(rc, TREXIO_SUCCESS)
    end if
  end if



! One-e RDM
! ---------

  rc = trexio_write_rdm_1e(f,one_e_dm_mo)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_rdm_1e_up(f,one_e_dm_mo_alpha_average)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_rdm_1e_dn(f,one_e_dm_mo_beta_average)
  call trexio_assert(rc, TREXIO_SUCCESS)


! Two-e RDM
! ---------

  if (export_rdm) then
    PROVIDE two_e_dm_mo
    print *, 'Two-e RDM'

    if (cholesky_rdm) then

      double precision :: xnorm0, x
      integer :: io
      integer(4) :: chol_index(3,BUFSIZE)

      xnorm0 = 0.d0
      x  = 0.d0
      io = 0
      do l=1,mo_num
        do k=1,mo_num
          do j=1,mo_num
            do i=1,mo_num
                if(i.eq.k.and.j.eq.l) then
                   xnorm0 = xnorm0 + two_e_dm_mo(i,j,k,l)
                   if(two_e_dm_mo(i,j,k,l).le.cholesky_tolerance) then
                      x = x + two_e_dm_mo(i,j,k,l)
                      io = io + 1
                   endif
                endif
              enddo
           enddo
        enddo
      enddo
      rank = mo_num*mo_num-io
      write(6,'(/,2X,''******    CHOLESKY DECOMPOSITION OF 2-RDM    ******'')')
      write(6,'(2X,''Estimated number of Chol vectors'',I5)') rank
      write(6,'(2X,''Estimated error in norm'',F12.4)') -x
     
      integer :: rank
      double precision, allocatable :: two_e_dm_mo_cholesky(:,:,:)
      allocate(two_e_dm_mo_cholesky(mo_num, mo_num, mo_num*mo_num-io))

      call pivoted_cholesky(two_e_dm_mo, rank, cholesky_tolerance, mo_num**2, two_e_dm_mo_cholesky)

      write(6,'(/,2X,''Cholesky decomposition of active 2-RDM completed. Tol: '',E10.2)')cholesky_tolerance
      write(6,'(2X,''Number of Cholesky vectors: '',I4)')rank
      write(6,'(2X,''Number of orbitals and the NChol/mo_num ratio:'',I4,F8.2)')mo_num,float(rank)/float(mo_num)
 

      rc = trexio_write_rdm_chol_num(f, rank)
      call trexio_assert(rc, TREXIO_SUCCESS)

      icount = 0_8
      offset = 0_8
      do k=1, rank
        do j=1,mo_num
          do i=1,mo_num
            integral = two_e_dm_mo_cholesky(i,j,k)
            if (integral == 0.d0) cycle
            icount += 1_8
            eri_buffer(icount) = integral
            chol_index(1,icount) = i
            chol_index(2,icount) = j
            chol_index(3,icount) = k
            if (icount == BUFSIZE) then
              rc = trexio_write_rdm_2e_cholesky(f, offset, icount, chol_index, eri_buffer)
              call trexio_assert(rc, TREXIO_SUCCESS)
              offset += icount
              icount = 0_8
            end if
          end do
        end do
      end do
   
      if (icount > 0_8) then
        rc = trexio_write_rdm_2e_cholesky(f, offset, icount, chol_index, eri_buffer)
        call trexio_assert(rc, TREXIO_SUCCESS)
      end if

    else

      icount = 0_8
      offset = 0_8
      do l=1,mo_num
        do k=1,mo_num
          do j=1,mo_num
            do i=1,mo_num
              integral = two_e_dm_mo(i,j,k,l)
              if (integral == 0.d0) cycle
              icount += 1_8
              eri_buffer(icount) = integral
              eri_index(1,icount) = i
              eri_index(2,icount) = j
              eri_index(3,icount) = k
              eri_index(4,icount) = l
              if (icount == BUFSIZE) then
                rc = trexio_write_rdm_2e(f, offset, icount, eri_index, eri_buffer)
                call trexio_assert(rc, TREXIO_SUCCESS)
                offset += icount
                icount = 0_8
              end if
            end do
          end do
        end do
      end do
   
      if (icount > 0_8) then
        rc = trexio_write_rdm_2e(f, offset, icount, eri_index, eri_buffer)
        call trexio_assert(rc, TREXIO_SUCCESS)
      end if
    end if

  end if


! ------------------------------------------------------------------------------

  ! Determinants
  ! ------------

    integer*8, allocatable :: det_buffer(:,:,:)
    double precision, allocatable :: coef_buffer(:,:)
    integer :: nint

!    rc = trexio_read_determinant_int64_num(f, nint)
!    call trexio_assert(rc, TREXIO_SUCCESS)
    nint = N_int
    if (nint /= N_int) then
       stop 'Problem with N_int'
    endif
    allocate ( det_buffer(nint, 2, BUFSIZE), coef_buffer(BUFSIZE, n_states) )

    icount = 0_8
    offset = 0_8
    rc = trexio_write_state_num(f, n_states)
    call trexio_assert(rc, TREXIO_SUCCESS)

    rc = trexio_set_state (f, 0)
    call trexio_assert(rc, TREXIO_SUCCESS)
    do k=1,n_det
       icount += 1_8
       det_buffer(1:nint, 1:2, icount) = psi_det(1:N_int, 1:2, k)
       coef_buffer(icount,1:N_states) = psi_coef(k,1:N_states)
       if (icount == BUFSIZE) then
         call trexio_assert(rc, TREXIO_SUCCESS)
         rc = trexio_write_determinant_list(f, offset, icount, det_buffer)
         call trexio_assert(rc, TREXIO_SUCCESS)
         do i=1,N_states
           rc = trexio_set_state (f, i-1)
           call trexio_assert(rc, TREXIO_SUCCESS)
           rc = trexio_write_determinant_coefficient(f, offset, icount, coef_buffer(1,i))
         end do
         rc = trexio_set_state (f, 0)
         offset += icount
         icount = 0_8
       end if
    end do

  if (icount >= 0_8) then
       call trexio_assert(rc, TREXIO_SUCCESS)
       rc = trexio_write_determinant_list(f, offset, icount, det_buffer)
       call trexio_assert(rc, TREXIO_SUCCESS)
       do i=1,N_states
         rc = trexio_set_state (f, i-1)
         call trexio_assert(rc, TREXIO_SUCCESS)
         rc = trexio_write_determinant_coefficient(f, offset, icount, coef_buffer(1,i))
       end do
       rc = trexio_set_state (f, 0)
  end if

  deallocate ( det_buffer, coef_buffer )

  rc = trexio_close(f)
  call trexio_assert(rc, TREXIO_SUCCESS)

end


! -*- mode: f90 -*-
