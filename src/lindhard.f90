module lindhard_module
   use sparse
   use wmpi
   use para
   implicit none
contains

   subroutine compute_lindhard()
      implicit none

      integer :: nq1, nq2, nq3, nk1, nk2, nk3
      integer :: nq_total, nk_total
      integer :: iq, ik, iqx, iqy, iqz, ikx, iky, ikz
      integer :: neval, nvecs, ndimq, nnz, nnzmax, ierr
      integer :: iomega, ib, jb, iorb
      integer :: nomega, matrix_switch
      integer :: unit_id
      real(dp) :: dk_weight
      real(dp) :: Beta_fake
      real(dp) :: broad
      real(dp) :: k(3), q(3), kq(3)
      real(dp) :: occ_k, occ_kq
      real(dp) :: freq
      real(dp) :: phase
      real(dp) :: Hartree2eV
      real(dp) :: q_cart(3)
      real(dp), allocatable :: omega(:)
      real(dp), allocatable :: eig_k(:), eig_kq(:)
      real(dp), allocatable :: work_eigs(:)
      complex(dp), allocatable :: chi_local(:, :)
      complex(dp), allocatable :: chi_total(:, :)
      complex(dp), allocatable :: acoo_k(:), acoo_q(:)
      integer, allocatable :: icoo_k(:), jcoo_k(:)
      integer, allocatable :: icoo_q(:), jcoo_q(:)
      complex(dp), allocatable :: Hamk_dense(:, :)
      complex(dp), allocatable :: Hamkq_dense(:, :)
      real(dp), allocatable :: eig_dense(:)
      complex(dp), allocatable :: eigvec_work(:, :)
      complex(dp), allocatable :: eigvec_k(:, :)
      complex(dp), allocatable :: eigvec_kq(:, :)
      complex(dp) :: denominator
      complex(dp) :: overlap_amp
      complex(dp) :: phase_factor
      real(dp) :: overlap_weight
      logical :: ritzvec
      logical :: is_static
      logical :: use_sparse
      real(dp), external :: fermi
      complex(dp) :: sigma
      real(dp), parameter :: twopi = 6.283185307179586477_dp
      character(len=32), parameter :: fmt_header = '(a, 3f16.6)'
      character(len=32), parameter :: fmt_line = '(f16.8, 2f20.10)'

      nq1 = Lindhard_Nq1
      nq2 = Lindhard_Nq2
      nq3 = Lindhard_Nq3
      nk1 = Lindhard_Nk1
      nk2 = Lindhard_Nk2
      nk3 = Lindhard_Nk3

      nq_total = nq1*nq2*nq3
      nk_total = nk1*nk2*nk3

      if (nq_total <= 0 .or. nk_total <= 0) then
         if (cpuid == 0) then
            write(stdout, '(a)') '>> Lindhard calculation skipped because the k/q mesh is empty.'
         end if
         return
      end if

      ndimq = Num_wann
      neval = NumSelectedEigenVals
      if (neval <= 0) neval = ndimq
      if (ndimq > 2 .and. neval > ndimq - 2) neval = ndimq - 2
      if (neval < 1) neval = min(ndimq, max(1, ndimq - 2))

      nvecs = int(1.5_dp*neval)
      if (nvecs < 50) nvecs = 50
      if (nvecs > ndimq) nvecs = ndimq

      nnz = splen
      nnzmax = splen + ndimq

      use_sparse = Is_Sparse .or. Is_Sparse_Hr

      nomega = Lindhard_omega_num
      if (nomega < 1) nomega = 1
      is_static = (trim(Lindhard_mode) == 'STATIC')
      if (is_static) nomega = 1
      Hartree2eV = 1.0_dp / eV2Hartree

      matrix_switch = 1
      select case (trim(Lindhard_matrix))
      case ('CMA')
         matrix_switch = 0
      case ('FULL')
         if (.not. allocated(Origin_cell%wannier_centers_direct)) then
            if (cpuid == 0) write(stdout, '(a)') 'WARNING: Wannier centers not available; falling back to DIAG Lindhard matrix elements.'
            matrix_switch = 1
         else
            if (size(Origin_cell%wannier_centers_direct, 2) < ndimq) then
               if (cpuid == 0) write(stdout, '(a)') 'WARNING: Incomplete Wannier centers; falling back to DIAG Lindhard matrix elements.'
               matrix_switch = 1
            else
            matrix_switch = 2
            end if
         end if
      case default
         matrix_switch = 1
      end select

      allocate(omega(nomega))
      if (is_static) then
         omega(1) = 0.0_dp
      else
         if (nomega <= 1) then
            omega(1) = Lindhard_omega_min
         else
            do iomega = 1, nomega
            omega(iomega) = Lindhard_omega_min + (Lindhard_omega_max - Lindhard_omega_min) * &
               (dble(iomega - 1) / dble(max(1, nomega - 1)))
            end do
         end if
      end if

      allocate(eig_k(neval))
      allocate(eig_kq(neval))
      allocate(work_eigs(neval))
      allocate(eigvec_work(ndimq, nvecs))
      allocate(eigvec_k(ndimq, neval))
      allocate(eigvec_kq(ndimq, neval))
      if (use_sparse) then
         allocate(acoo_k(nnzmax))
         allocate(icoo_k(nnzmax))
         allocate(jcoo_k(nnzmax))
         allocate(acoo_q(nnzmax))
         allocate(icoo_q(nnzmax))
         allocate(jcoo_q(nnzmax))
      else
         allocate(Hamk_dense(ndimq, ndimq))
         allocate(Hamkq_dense(ndimq, ndimq))
         allocate(eig_dense(ndimq))
      end if
      allocate(chi_local(nomega, nq_total))
      allocate(chi_total(nomega, nq_total))

      chi_local = (0.0_dp, 0.0_dp)
      chi_total = (0.0_dp, 0.0_dp)

      dk_weight = 1.0_dp / dble(nk_total)
      Beta_fake = Beta
      if (Beta_fake <= 0.0_dp) Beta_fake = 1.0d12
      broad = Lindhard_broadening

      ritzvec = .true.
      sigma = cmplx(iso_energy, 0.0_dp, kind=dp)

      do iq = 1 + cpuid, nq_total, num_cpu
         iqx = (iq - 1) / (nq2*nq3) + 1
         iqy = ((iq - 1 - (iqx - 1)*nq2*nq3) / nq3) + 1
         iqz = iq - (iqy - 1)*nq3 - (iqx - 1)*nq2*nq3

         q = Lindhard_q_start + K3D_vec1_cube * dble(iqx - 1) / dble(max(1, nq1)) &
            + K3D_vec2_cube * dble(iqy - 1) / dble(max(1, nq2)) &
            + K3D_vec3_cube * dble(iqz - 1) / dble(max(1, nq3))

         do ik = 1, nk_total
            ikx = (ik - 1) / (nk2*nk3) + 1
            iky = ((ik - 1 - (ikx - 1)*nk2*nk3) / nk3) + 1
            ikz = ik - (iky - 1)*nk3 - (ikx - 1)*nk2*nk3

            k = K3D_start_cube + K3D_vec1_cube * dble(ikx - 1) / dble(max(1, nk1)) &
               + K3D_vec2_cube * dble(iky - 1) / dble(max(1, nk2)) &
               + K3D_vec3_cube * dble(ikz - 1) / dble(max(1, nk3))

            kq = k + q
            kq(1) = modulo(kq(1), 1.0_dp)
            kq(2) = modulo(kq(2), 1.0_dp)
            kq(3) = modulo(kq(3), 1.0_dp)

            if (use_sparse) then
               nnz = splen
               call ham_bulk_coo_sparsehr(k, acoo_k, icoo_k, jcoo_k)
               call arpack_sparse_coo_eigs(ndimq, nnzmax, nnz, acoo_k, jcoo_k, icoo_k, neval, nvecs, &
                  work_eigs, sigma, eigvec_work, ritzvec)
               eig_k(:) = work_eigs
               eigvec_k(:, :) = eigvec_work(:, 1:neval)

               nnz = splen
               call ham_bulk_coo_sparsehr(kq, acoo_q, icoo_q, jcoo_q)
               call arpack_sparse_coo_eigs(ndimq, nnzmax, nnz, acoo_q, jcoo_q, icoo_q, neval, nvecs, &
                  work_eigs, sigma, eigvec_work, ritzvec)
               eig_kq(:) = work_eigs
               eigvec_kq(:, :) = eigvec_work(:, 1:neval)
            else
               Hamk_dense = (0.0_dp, 0.0_dp)
               Hamkq_dense = (0.0_dp, 0.0_dp)
               call ham_bulk_atomicgauge(k, Hamk_dense)
               call eigensystem_c('V', 'U', ndimq, Hamk_dense, eig_dense)
               eig_k(:) = eig_dense(1:neval)
               eigvec_k(:, :) = Hamk_dense(:, 1:neval)

               call ham_bulk_atomicgauge(kq, Hamkq_dense)
               call eigensystem_c('V', 'U', ndimq, Hamkq_dense, eig_dense)
               eig_kq(:) = eig_dense(1:neval)
               eigvec_kq(:, :) = Hamkq_dense(:, 1:neval)
            end if

            do ib = 1, neval
               occ_k = fermi(eig_k(ib) - iso_energy, Beta_fake)
               do jb = 1, neval
                  occ_kq = fermi(eig_kq(jb) - iso_energy, Beta_fake)
                  if (abs(occ_k - occ_kq) < 1.0d-12) cycle
                  select case (matrix_switch)
                  case (0)
                     overlap_weight = 1.0_dp
                  case (1)
                     overlap_amp = dot_product(conjg(eigvec_k(:, ib)), eigvec_kq(:, jb))
                     overlap_weight = abs(overlap_amp)**2
                  case (2)
                     overlap_amp = (0.0_dp, 0.0_dp)
                     do iorb = 1, ndimq
                        phase = twopi * (q(1)*Origin_cell%wannier_centers_direct(1, iorb) + &
                                         q(2)*Origin_cell%wannier_centers_direct(2, iorb) + &
                                         q(3)*Origin_cell%wannier_centers_direct(3, iorb))
                        phase_factor = cmplx(cos(phase), -sin(phase), kind=dp)
                        overlap_amp = overlap_amp + conjg(eigvec_k(iorb, ib))*eigvec_kq(iorb, jb)*phase_factor
                     end do
                     overlap_weight = abs(overlap_amp)**2
                  end select
                  if (overlap_weight <= 0.0_dp) cycle
                  do iomega = 1, nomega
                     freq = omega(iomega)
                     denominator = (freq + zi*broad) + (eig_k(ib) - eig_kq(jb))
                     if (abs(denominator) <= 1.0d-14) cycle
                     chi_local(iomega, iq) = chi_local(iomega, iq) + (occ_k - occ_kq) * overlap_weight / denominator
                  end do
               end do
            end do
         end do

         chi_local(:, iq) = chi_local(:, iq) * dk_weight
      end do

#if defined (MPI)
      call mpi_allreduce(chi_local, chi_total, size(chi_local), mpi_dc, mpi_sum, mpi_cmw, ierr)
#else
      chi_total = chi_local
#endif

      if (cpuid == 0) then
         unit_id = outfileindex + 1
         outfileindex = outfileindex + 1
         open(unit=unit_id, file=trim(Lindhard_output))
         write(unit_id, '(a)') '# Lindhard susceptibility chi(q, omega)'
         write(unit_id, '(a)') '# omega in eV, chi in 1/eV'
         write(unit_id, '(a, a)') '# mode = ', trim(Lindhard_mode)
         write(unit_id, '(a, a)') '# matrix = ', trim(Lindhard_matrix)
         if (Lindhard_T > 0.0_dp) then
            write(unit_id, '(a, f16.6, 1x, a)') '# temperature = ', Lindhard_T, trim(Lindhard_T_unit)
         end if
         write(unit_id, '(a, f16.8)') '# beta (Hartree^-1) = ', Beta_fake
         write(unit_id, '(a, f16.8)') '# broadening (eV) = ', Lindhard_broadening*Hartree2eV

         do iq = 1, nq_total
            iqx = (iq - 1) / (nq2*nq3) + 1
            iqy = ((iq - 1 - (iqx - 1)*nq2*nq3) / nq3) + 1
            iqz = iq - (iqy - 1)*nq3 - (iqx - 1)*nq2*nq3

            q = Lindhard_q_start + K3D_vec1_cube * dble(iqx - 1) / dble(max(1, nq1)) &
               + K3D_vec2_cube * dble(iqy - 1) / dble(max(1, nq2)) &
               + K3D_vec3_cube * dble(iqz - 1) / dble(max(1, nq3))
            write(unit_id, fmt_header) '# q (fractional) =', q
            q_cart = q(1)*Origin_cell%Kua + q(2)*Origin_cell%Kub + q(3)*Origin_cell%Kuc
            write(unit_id, fmt_header) '# q (cartesian 1/Angstrom) =', q_cart*Angstrom2atomic
            do iomega = 1, nomega
               write(unit_id, fmt_line) omega(iomega)/eV2Hartree, &
                  real(chi_total(iomega, iq))*eV2Hartree, aimag(chi_total(iomega, iq))*eV2Hartree
            end do
            write(unit_id, '(a)') ' '
         end do
         close(unit_id)
      end if

      deallocate(omega)
      deallocate(eig_k)
      deallocate(eig_kq)
      deallocate(work_eigs)
      deallocate(eigvec_work)
      deallocate(eigvec_k)
      deallocate(eigvec_kq)
      if (use_sparse) then
         deallocate(acoo_k)
         deallocate(icoo_k)
         deallocate(jcoo_k)
         deallocate(acoo_q)
         deallocate(icoo_q)
         deallocate(jcoo_q)
      else
         deallocate(Hamk_dense)
         deallocate(Hamkq_dense)
         deallocate(eig_dense)
      end if
      deallocate(chi_local)
      deallocate(chi_total)

   end subroutine compute_lindhard

end module lindhard_module
