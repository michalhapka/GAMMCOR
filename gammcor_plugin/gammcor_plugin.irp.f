program gammcor_plugin
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .True.
  TOUCH read_wf 
  call run
end

subroutine run
  implicit none

  print *, 'nuclear_repulsion: ', nuclear_repulsion
  call print_energy_components_local 

end



subroutine print_energy_components_local
  implicit none
  BEGIN_DOC
! Prints the different components of the energy.
  END_DOC
  integer, save :: ifirst = 0
  double precision :: Vee, Ven, Vnn, Vecp, T, f
  integer  :: i,j,k

  Vnn = nuclear_repulsion

  print *, 'Energy components'
  print *, '================='
  print *, ''
  do k=1,N_states

    Ven  = 0.d0
    Vecp = 0.d0
    T    = 0.d0

    do j=1,mo_num
      do i=1,mo_num
        f = one_e_dm_mo_alpha(i,j,k) + one_e_dm_mo_beta(i,j,k)
        Ven  = Ven  + f * mo_integrals_n_e(i,j)
        Vecp = Vecp + f * mo_pseudo_integrals(i,j)
        T    = T    + f * mo_kinetic_integrals(i,j)
      enddo
    enddo
    Vee = psi_energy(k) - Ven - Vecp - T
    
    if (ifirst == 0) then
      ifirst = 1
      print *, 'Vnn  : Nucleus-Nucleus   potential energy'
      print *, 'Ven  : Electron-Nucleus  potential energy'
      print *, 'Vee  : Electron-Electron potential energy'
      print *, 'Vecp : Potential energy of the pseudo-potentials'
      print *, 'T    : Electronic kinetic energy'
      print *, ''
    endif

    print *, 'State ', k
    print *, '---------'
    print *, ''
    print *, 'Vnn  = ', Vnn
    print *, 'One-e  ', Ven+T
    print *, 'Vee  = ', Vee
    print *, 'Eele = ', Ven + T + Vee
    print *, 'ETot = ', Vnn + Ven + T + Vee
    print *, ''
  enddo

  print *, ''

end
