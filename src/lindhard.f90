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
      integer :: iomega, ib, jb
      real(dp) :: dk_weight
      real(dp) :: Beta_fake
      real(dp) :: broad
      real(dp) :: k(3), q(3), kq(3)
      real(dp) :: occ_k, occ_kq
      real(dp) :: freq
      real(dp), allocatable :: omega(:)
      real(dp), allocatable :: eig_k(:), eig_kq(:)
      real(dp), allocatable :: work_eigs(:)
      complex(dp), allocatable :: chi_local(:, :)
      complex(dp), allocatable :: chi_total(:, :)
      complex(dp), allocatable :: acoo_k(:), acoo_q(:)
      integer, allocatable :: icoo_k(:), jcoo_k(:)
      integer, allocatable :: icoo_q(:), jcoo_q(:)
      complex(dp), allocatable :: zeigv(:, :)
      complex(dp) :: denom
      logical :: ritzvec
      real(dp), external :: fermi

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
      if (ndimq > 2) then
         if (neval > ndimq - 2) neval = ndimq - 2
      else
         neval = ndimq
      end if
      if (neval < 1) neval = ndimq

      nvecs = int(1.5_dp*neval)
      if (nvecs < 50) nvecs = 50
      if (nvecs > ndimq) nvecs = ndimq

      nnz = splen
      nnzmax = splen + ndimq

      allocate(omega(Lindhard_omega_num))
      if (Lindhard_omega_num <= 1) then
         omega(1) = Lindhard_omega_min
      else
         do iomega = 1, Lindhard_omega_num
            omega(iomega) = Lindhard_omega_min + (Lindhard_omega_max - Lindhard_omega_min) * &
               (dble(iomega - 1) / dble(Lindhard_omega_num - 1))
         end do
      end if

      allocate(eig_k(neval))
      allocate(eig_kq(neval))
      allocate(work_eigs(ndimq))
      allocate(zeigv(ndimq, nvecs))
      allocate(acoo_k(nnzmax))
      allocate(icoo_k(nnzmax))
      allocate(jcoo_k(nnzmax))
      allocate(acoo_q(nnzmax))
      allocate(icoo_q(nnzmax))
      allocate(jcoo_q(nnzmax))
      allocate(chi_local(Lindhard_omega_num, nq_total))
      allocate(chi_total(Lindhard_omega_num, nq_total))

      chi_local = (0.0_dp, 0.0_dp)
      chi_total = (0.0_dp, 0.0_dp)

      dk_weight = 1.0_dp / dble(nk_total)
      Beta_fake = Beta
      if (Beta_fake <= 0.0_dp) Beta_fake = 1.0d12
      broad = Lindhard_broadening

      ritzvec = .false.

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

            call ham_bulk_coo_sparsehr(k, acoo_k, icoo_k, jcoo_k)
            call arpack_sparse_coo_eigs(ndimq, nnzmax, nnz, acoo_k, jcoo_k, icoo_k, neval, nvecs, &
               work_eigs, (0.0_dp, 0.0_dp), zeigv, ritzvec)
            eig_k(:) = work_eigs(1:neval)

            call ham_bulk_coo_sparsehr(kq, acoo_q, icoo_q, jcoo_q)
            call arpack_sparse_coo_eigs(ndimq, nnzmax, nnz, acoo_q, jcoo_q, icoo_q, neval, nvecs, &
               work_eigs, (0.0_dp, 0.0_dp), zeigv, ritzvec)
            eig_kq(:) = work_eigs(1:neval)

            do ib = 1, neval
               occ_k = fermi(eig_k(ib) - iso_energy, Beta_fake)
               do jb = 1, neval
                  occ_kq = fermi(eig_kq(jb) - iso_energy, Beta_fake)
                  if (abs(occ_k - occ_kq) < 1.0d-12) cycle
                  do iomega = 1, Lindhard_omega_num
                     freq = omega(iomega)
                     denom = (freq + zi*broad) + (eig_k(ib) - eig_kq(jb))
                     chi_local(iomega, iq) = chi_local(iomega, iq) + (occ_k - occ_kq) / denom
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
         character(len=*), parameter :: fmt_header = '(a, 3f16.6)'
         character(len=*), parameter :: fmt_line = '(f16.8, 2f20.10)'
         real(dp) :: q_cart(3)
         integer :: unit_id

         unit_id = outfileindex + 1
         outfileindex = outfileindex + 1
         open(unit=unit_id, file=trim(Lindhard_output))
         write(unit_id, '(a)') '# Lindhard susceptibility chi(q, omega)'
         write(unit_id, '(a)') '# omega in eV, chi in 1/eV'
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
            do iomega = 1, Lindhard_omega_num
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
      deallocate(zeigv)
      deallocate(acoo_k)
      deallocate(icoo_k)
      deallocate(jcoo_k)
      deallocate(acoo_q)
      deallocate(icoo_q)
      deallocate(jcoo_q)
      deallocate(chi_local)
      deallocate(chi_total)

   end subroutine compute_lindhard

end module lindhard_module
