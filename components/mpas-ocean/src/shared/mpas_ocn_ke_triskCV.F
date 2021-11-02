! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the
! LICENSE file
! distributed with this code, or at
! http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_ke_triskCV
!
!> \brief MPAS ocean driver for the kinetic energy in TRiSK-CV
!> \author Sara Calandrini
!> \date   17 September 2021
!> \details
!>  This module contains the routines for computing the kinetic 
!>  energy in TRiSK-CV. 
!
!-----------------------------------------------------------------------

module ocn_ke_triskCV

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_stream_manager
   use mpas_io_units
   use mpas_dmpar

   use ocn_constants
   use ocn_config 
   use ocn_thick_ale
   use ocn_diagnostics_variables
   use ocn_mesh

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: controlVolumeEdge, &
             buildInvMatrixBtB, &
             Wlsqr, &
             inv_3x3, &
             sphericalTriangleArea, &
             sphere_distance

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------  

!***********************************************************************

contains

   subroutine controlVolumeEdge(domain, meshPool)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! controlVolumeEdge computes the areas associated to every edge in the 
   ! TRiSK++ discretization 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (domain_type), intent(inout) :: domain
      integer :: vertex1, vertex2, cell1, cell2, iEdge, jEdge, i, k
      real (kind=RKIND) :: a, b 
      real (kind=RKIND), dimension(:), pointer :: areaCVEdge
      real (kind=RKIND), dimension(:,:), pointer :: weightsEdgeToCell

      call mpas_pool_get_array(meshPool, 'areaCVEdge', areaCVEdge)
      call mpas_pool_get_array(meshPool, 'weightsEdgeToCell', weightsEdgeToCell)


      areaCVEdge(:) = 0.0_RKIND
      weightsEdgeToCell(:,:) = 0.0_RKIND
      do iEdge = 1, nEdgesAll !nEdgesArray(3)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)
         a = sphericalTriangleArea(latVertex(vertex1), lonVertex(vertex1), latVertex(vertex2), lonVertex(vertex2), latCell(cell1), lonCell(cell1))
         !if (latCell(cell2)==0.0 .AND. lonCell(cell2)==0.0) then
         if (cell2 > nCellsAll) then 
            b = 0.0_RKIND
         else
            b = sphericalTriangleArea(latVertex(vertex1), lonVertex(vertex1), latVertex(vertex2), lonVertex(vertex2), latCell(cell2), lonCell(cell2))
         end if 
         areaCVEdge(iEdge) = a+b
         do i = 1, nEdgesOnCell(cell1)
            jEdge = edgesOnCell(i,cell1)
            if (iEdge == jEdge) then
               weightsEdgeToCell(i,cell1) = a
            end if
         end do
         do i = 1, nEdgesOnCell(cell2)
            jEdge = edgesOnCell(i,cell2)
            if (iEdge == jEdge) then
               weightsEdgeToCell(i,cell2) = b
            end if
         end do
      end do
      call mpas_dmpar_field_halo_exch(domain, 'weightsEdgeToCell')
      call mpas_dmpar_field_halo_exch(domain, 'areaCVEdge')

      !do iEdge = 1, nEdgesAll
      !   print*, iEdge, areaCVEdge(iEdge) 
      !end do

   end subroutine controlVolumeEdge

   subroutine buildInvMatrixBtB(domain, meshPool)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! buildMatrixBtB computes the least square matrix B^{T}B needed for the
   ! computation of the kinetic energy in the TRiSK++ discretization 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (domain_type), intent(inout) :: domain
      integer :: iEdge, jEdge, i, j, k
      integer, dimension(2) :: iEdgeCells, jEdgeCells
      real (kind=RKIND) :: re_x, re_y, re_z, ne_x, ne_y, ne_z, norm, xTangent, yTangent, zTangent, detinv
      real (kind=RKIND), dimension(:), pointer :: areaCVEdge   
      real (kind=RKIND), dimension(:,:,:), pointer :: invBtB
      real (kind=RKIND), dimension(:,:), allocatable :: BtB

      call mpas_pool_get_array(meshPool, 'areaCVEdge', areaCVEdge)
      !call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
      !call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'invBtB', invBtB)

      allocate(BtB(vertexDegree, vertexDegree))

      invBtB(:,:,:) = 0.0_RKIND
      !print*, nEdgesHalo

      !do iEdge = 1, nEdgesAll
      !   print*, iEdge, areaCVEdge(iEdge)
      !end do

      do iEdge = 1, nEdgesHalo(3)
         BtB(:,:) = 0.0_RKIND
         iEdgeCells(1) = cellsOnEdge(1,iEdge)
         iEdgeCells(2) = cellsOnEdge(2,iEdge)
         norm = sqrt(xEdge(iEdge)**2.0 + yEdge(iEdge)**2.0 + zEdge(iEdge)**2.0)
         re_x = xEdge(iEdge)/norm
         re_y = yEdge(iEdge)/norm
         re_z = zEdge(iEdge)/norm
         !print*, re_x*ne_x + re_y*ne_y + re_z*ne_z !e-3
         do j = 1, nEdgesOnEdge(iEdge)
            !print*, nEdgesOnEdge(iEdge)
            jEdge = edgesOnEdge(j,iEdge)
            jEdgeCells(1) = cellsOnEdge(1,jEdge)
            jEdgeCells(2) = cellsOnEdge(2,jEdge)
            ne_x = (xCell(jEdgeCells(1)) - xCell(jEdgeCells(2)))
            ne_y = (yCell(jEdgeCells(1)) - yCell(jEdgeCells(2)))
            ne_z = (zCell(jEdgeCells(1)) - zCell(jEdgeCells(2)))
            norm = sqrt(ne_x**2.0 + ne_y**2.0 + ne_z**2.0)
            ne_x = ne_x/norm
            ne_y = ne_y/norm
            ne_z = ne_z/norm
            BtB(1,1) = BtB(1,1) + areaCVEdge(jEdge) * ne_x * areaCVEdge(jEdge) * ne_x
            BtB(2,2) = BtB(2,2) + areaCVEdge(jEdge) * ne_y * areaCVEdge(jEdge) * ne_y
            BtB(3,3) = BtB(3,3) + areaCVEdge(jEdge) * ne_z * areaCVEdge(jEdge) * ne_z
            BtB(1,2) = BtB(1,2) + areaCVEdge(jEdge) * ne_x * areaCVEdge(jEdge) * ne_y
            BtB(1,3) = BtB(1,3) + areaCVEdge(jEdge) * ne_x * areaCVEdge(jEdge) * ne_z
            BtB(2,3) = BtB(2,3) + areaCVEdge(jEdge) * ne_y * areaCVEdge(jEdge) * ne_z
         end do
         BtB(1,1) = BtB(1,1) + re_x*re_x
         BtB(2,2) = BtB(2,2) + re_y*re_y
         BtB(3,3) = BtB(3,3) + re_z*re_z
         BtB(1,2) = BtB(1,2) + re_x*re_y
         BtB(1,3) = BtB(1,3) + re_x*re_z
         BtB(2,1) = BtB(1,2)
         BtB(2,3) = BtB(2,3) + re_y*re_z
         BtB(3,1) = BtB(1,3)
         BtB(3,2) = BtB(2,3)

         !Calculate the inverse determinant of the matrix

         !print*, iEdge, BtB(1,1), BtB(2,2), BtB(3,3), BtB(1,2), BtB(1,3), BtB(2,1), &
         !        BtB(2,3), BtB(3,1), BtB(3,2)

         !print*, iEdge, nEdgesOnEdge(iEdge), iEdgeCells(1), iEdgeCells(2), (BtB(1,1)*BtB(2,2)*BtB(3,3) - BtB(1,1)*BtB(2,3)*BtB(3,2)&
         !       - BtB(1,2)*BtB(2,1)*BtB(3,3) + BtB(1,2)*BtB(2,3)*BtB(3,1)&
         !       + BtB(1,3)*BtB(2,1)*BtB(3,2) - BtB(1,3)*BtB(2,2)*BtB(3,1))

         if (iEdgeCells(2) > nCellsAll) then 
            detinv = 0.0_RKIND
         else
            detinv = 1.0_RKIND/(BtB(1,1)*BtB(2,2)*BtB(3,3) - BtB(1,1)*BtB(2,3)*BtB(3,2)&
                   - BtB(1,2)*BtB(2,1)*BtB(3,3) + BtB(1,2)*BtB(2,3)*BtB(3,1)&
                   + BtB(1,3)*BtB(2,1)*BtB(3,2) - BtB(1,3)*BtB(2,2)*BtB(3,1))
         end if
         !print*, iEdge, detinv

         ! Calculate the inverse of the matrix
         invBtB(1,1,iEdge) = +detinv * (BtB(2,2)*BtB(3,3) - BtB(2,3)*BtB(3,2))
         invBtB(2,1,iEdge) = -detinv * (BtB(2,1)*BtB(3,3) - BtB(2,3)*BtB(3,1))
         invBtB(3,1,iEdge) = +detinv * (BtB(2,1)*BtB(3,2) - BtB(2,2)*BtB(3,1))
         invBtB(1,2,iEdge) = -detinv * (BtB(1,2)*BtB(3,3) - BtB(1,3)*BtB(3,2))
         invBtB(2,2,iEdge) = +detinv * (BtB(1,1)*BtB(3,3) - BtB(1,3)*BtB(3,1))
         invBtB(3,2,iEdge) = -detinv * (BtB(1,1)*BtB(3,2) - BtB(1,2)*BtB(3,1))
         invBtB(1,3,iEdge) = +detinv * (BtB(1,2)*BtB(2,3) - BtB(1,3)*BtB(2,2))
         invBtB(2,3,iEdge) = -detinv * (BtB(1,1)*BtB(2,3) - BtB(1,3)*BtB(2,1))
         invBtB(3,3,iEdge) = +detinv * (BtB(1,1)*BtB(2,2) - BtB(1,2)*BtB(2,1))
      end do
      call mpas_dmpar_field_halo_exch(domain, 'invBtB')

      deallocate(BtB)

   end subroutine buildInvMatrixBtB

   subroutine inv_3x3(amat,adim,ainv,vdim, adet)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! inv3x3 computes the inverse of a 3x3 matrix 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       implicit none

   !------------------------------------------- arguments !
       integer, intent( in) :: adim
       real (kind=RKIND), intent( in) :: amat(adim,*)
       integer, intent( in) :: vdim
       real (kind=RKIND), intent(out) :: ainv(vdim,*)
       real (kind=RKIND), intent(out) :: adet

   !------------------------------------------- variables !
       real (kind=RKIND) :: &
       aa2233,aa2332,aa2133,aa2331,aa2132, &
       aa2231,aa1233,aa1332,aa1223,aa1322, &
       aa1133,aa1331,aa1123,aa1321,aa1132, &
       aa1231,aa1122,aa1221

   !------------------------------------------- form A^-1 !
       aa2233 = amat(2,2) * amat(3,3)
       aa2332 = amat(2,3) * amat(3,2)
       aa2133 = amat(2,1) * amat(3,3)
       aa2331 = amat(2,3) * amat(3,1)
       aa2132 = amat(2,1) * amat(3,2)
       aa2231 = amat(2,2) * amat(3,1)

       adet =  &
       amat(1,1) *  (aa2233 - aa2332) - &
           amat(1,2) *  (aa2133 - aa2331) + &
           amat(1,3) *  (aa2132 - aa2231)

       aa1233 = amat(1,2) * amat(3,3)
       aa1332 = amat(1,3) * amat(3,2)
       aa1223 = amat(1,2) * amat(2,3)
       aa1322 = amat(1,3) * amat(2,2)
       aa1133 = amat(1,1) * amat(3,3)
       aa1331 = amat(1,3) * amat(3,1)
       aa1123 = amat(1,1) * amat(2,3)
       aa1321 = amat(1,3) * amat(2,1)
       aa1132 = amat(1,1) * amat(3,2)
       aa1231 = amat(1,2) * amat(3,1)
       aa1122 = amat(1,1) * amat(2,2)
       aa1221 = amat(1,2) * amat(2,1)

       ainv(1,1) =  (aa2233 - aa2332)
       ainv(1,2) = -(aa1233 - aa1332)
       ainv(1,3) =  (aa1223 - aa1322)

       ainv(2,1) = -(aa2133 - aa2331)
       ainv(2,2) =  (aa1133 - aa1331)
       ainv(2,3) = -(aa1123 - aa1321)

       ainv(3,1) =  (aa2132 - aa2231)
       ainv(3,2) = -(aa1132 - aa1231)
       ainv(3,3) =  (aa1122 - aa1221)

       return

   end subroutine inv_3x3

   subroutine Wlsqr(domain, meshPool)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! buildMatrixBtB computes the least square weightsOnEdge needed for the
   ! computation of the kinetic energy in the TRiSK++ discretization 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (domain_type), intent(inout) :: domain
      integer :: iEdge, jEdge, iEoE, i, j, k, cell1
      integer, dimension(2) :: iEdgeCells, jEdgeCells
      real (kind=RKIND) :: re_x, re_y, re_z, ne_x, ne_y, ne_z, norm, xTangent, yTangent, zTangent, normTangent, detinv   
      real (kind=RKIND), dimension(:), pointer :: areaCVEdge
      real (kind=RKIND), dimension(:,:,:), pointer :: invBtB

      real (kind=RKIND), pointer :: sphere_radius
      real (kind=RKIND) :: tmat(25,3), matt(3,25), rmat(3,3), rinv(3,3), &
                           rdet, xmul, ymul, zmul
      real (kind=RKIND), dimension(:), allocatable :: elen
      real (kind=RKIND), dimension(:,:), allocatable :: enrm, edir, eprp

      call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)
      call mpas_pool_get_array(meshPool, 'areaCVEdge', areaCVEdge) 
      call mpas_pool_get_array(meshPool, 'invBtB', invBtB)
      !call mpas_pool_get_array(meshPool, 'weightsOnEdge_lsqr', weightsOnEdge_lsqr)
      !call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)

      allocate(enrm(nEdgesHalo(3),3))
      allocate(edir(nEdgesAll,3))
      allocate(eprp(nEdgesHalo(3),3))
      allocate(elen(nEdgesAll))

   !- build lsqr weights for each edge in the mesh
      do iEdge = 1, nEdgesHalo(3)
         tmat(:,:) = 0.0_RKIND
         matt(:,:) = 0.0_RKIND
         rmat(:,:) = 0.0_RKIND
         rinv(:,:) = 0.0_RKIND

      !- unit tangent (on spherical surface) to each edge (in E^3)
         eprp(iEdge,1) = xVertex(verticesOnEdge(2,iEdge)) - &
                    xVertex(verticesOnEdge(1,iEdge))
         eprp(iEdge,2) = yVertex(verticesOnEdge(2,iEdge)) - &
                    yVertex(verticesOnEdge(1,iEdge))
         eprp(iEdge,3) = zVertex(verticesOnEdge(2,iEdge)) - &
                    zVertex(verticesOnEdge(1,iEdge))

         norm = sqrt(eprp(iEdge,1)**2 + eprp(iEdge,2)**2 + eprp(iEdge,3)**2)

         eprp(iEdge,1) = eprp(iEdge,1) / norm
         eprp(iEdge,2) = eprp(iEdge,2) / norm
         eprp(iEdge,3) = eprp(iEdge,3) / norm

      !- unit normals (to spherical surface) at each edge
         enrm(iEdge,1) = xEdge(iEdge) / sphere_radius
         enrm(iEdge,2) = yEdge(iEdge) / sphere_radius
         enrm(iEdge,3) = zEdge(iEdge) / sphere_radius

      !- assemble lsqr matrix T: n_hat dot u = u_e
         tmat(1,1) = enrm(iEdge,1)
         tmat(1,2) = enrm(iEdge,2)
         tmat(1,3) = enrm(iEdge,3)

         do jEdge = 1, nEdgesOnEdge(iEdge)
            iEoE = edgesOnEdge(jEdge,iEdge)

         !- unit normals (on spherical surface) to each edge (in E^3)
            if (cellsOnEdge(2,iEoE)>nCellsAll) then
               cell1 = cellsOnEdge(1,iEoE)
               !do i = 1, nEdgesOnCell(cell1)
               !   if(iEoE == edgesOnCell(i,cell1)) then  
               !      edir(iEoE,1) = edgeSignOnCell(i,cell1)*(xEdge(iEoE) - xCell(cell1))
               !      edir(iEoE,2) = edgeSignOnCell(i,cell1)*(yEdge(iEoE) - yCell(cell1))
               !      edir(iEoE,3) = edgeSignOnCell(i,cell1)*(zEdge(iEoE) - zCell(cell1))
               !   end if
               !end do
               edir(iEoE,1) = xEdge(iEoE) - xCell(cell1)
               edir(iEoE,2) = yEdge(iEoE) - yCell(cell1)
               edir(iEoE,3) = zEdge(iEoE) - zCell(cell1)
            else
               edir(iEoE,1) = xCell(cellsOnEdge(2,iEoE)) - &
                            xCell(cellsOnEdge(1,iEoE))
               edir(iEoE,2) = yCell(cellsOnEdge(2,iEoE)) - &
                            yCell(cellsOnEdge(1,iEoE))
               edir(iEoE,3) = zCell(cellsOnEdge(2,iEoE)) - &
                            zCell(cellsOnEdge(1,iEoE))
            end if 

            elen(iEoE) = sqrt(edir(iEoE,1)**2 + edir(iEoE,2)**2 + edir(iEoE,3)**2)

           ! if(iEdge==2380) then
           !    print*, jEdge, cellsOnEdge(1,iEoE), cellsOnEdge(2,iEoE), xCell(cellsOnEdge(2,iEoE)), yCell(cellsOnEdge(2,iEoE)), zCell(cellsOnEdge(2,iEoE)), elen(iEoE)
           ! end if

            edir(iEoE,1) = edir(iEoE,1) / elen(iEoE)
            edir(iEoE,2) = edir(iEoE,2) / elen(iEoE)
            edir(iEoE,3) = edir(iEoE,3) / elen(iEoE)

            tmat(jEdge+1,1) = edir(iEoE,1)
            tmat(jEdge+1,2) = edir(iEoE,2)
            tmat(jEdge+1,3) = edir(iEoE,3)
         end do

      !- form lsqr matrices: (T'*T) * u = T' * u_e
         matt = transpose(tmat)

         rmat = matmul(matt,tmat)

      !- factorise matrices: r_inv = (T'*T) ^ -1
         call inv_3x3(rmat, 3, rinv, 3, rdet)

      !- assemble lsqr weight: t_hat dot r_inv * T'
         do jEdge = 1, nEdgesOnEdge(iEdge)
            xmul = ( &
               rinv(1,1) * matt(1,jEdge+1) + &
               rinv(1,2) * matt(2,jEdge+1) + &
               rinv(1,3) * matt(3,jEdge+1) &
            ) / rdet

            ymul = ( &
               rinv(2,1) * matt(1,jEdge+1) + &
               rinv(2,2) * matt(2,jEdge+1) + &
               rinv(2,3) * matt(3,jEdge+1) &
            ) / rdet

            zmul = ( &
               rinv(3,1) * matt(1,jEdge+1) + &
               rinv(3,2) * matt(2,jEdge+1) + &
               rinv(3,3) * matt(3,jEdge+1) &
            ) / rdet

            weightsOnEdge_lsqr(jEdge,iEdge) = &
               eprp(iEdge,1) * xmul + &
               eprp(iEdge,2) * ymul + &
               eprp(iEdge,3) * zmul

            if(iEdge==2380) then
            !   print*, jEdge, xmul, rinv(1,1), rinv(1,2), rinv(1,3), matt(1,jEdge+1), matt(2,jEdge+1), matt(3,jEdge+1), rdet 
               print*, jEdge, eprp(iEdge,1), eprp(iEdge,2), eprp(iEdge,3), xmul, ymul, zmul, weightsOnEdge_lsqr(jEdge,iEdge), weightsOnEdge(jEdge,iEdge)
            end if 

         end do
      end do

      call mpas_dmpar_field_halo_exch(domain, 'weightsOnEdge_lsqr')

      !do iEdge = 1, nEdgesAll
      !   do jEdge = 1, nEdgesOnEdge(iEdge)
      !      if (weightsOnEdge_lsqr(jEdge,iEdge) > 1.0_RKIND) then
      !      print*, iEdge, cellsOnEdge(1,iEdge), cellsOnEdge(2,iEdge), weightsOnEdge_lsqr(jEdge,iEdge), weightsOnEdge(jEdge,iEdge)
      !      end if
      !   end do
      !end do

      deallocate(eprp)
      deallocate(edir)
      deallocate(enrm)
      deallocate(elen) 

   end subroutine Wlsqr


   real function sphericalTriangleArea(laEdge, loEdge, laVertex, loVertex, laCell, loCell)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! sphericalTriangleArea uses the spherical analog of Heron's formula to
   ! compute the area of a triangle on the surface of a sphere
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), intent(in) :: laEdge, loEdge, laVertex, loVertex, laCell, loCell
      real (kind=RKIND) :: tanqe, s, a, b, c, one, radius
      one = 1.0
      radius = 6371220.0 

      a = sphere_distance(laCell, loCell, laVertex, loVertex, one)
      b = sphere_distance(laCell, loCell, laEdge, loEdge, one)
      c = sphere_distance(laEdge, loEdge, laVertex, loVertex, one)
      s = 0.5*(a+b+c)
      !print*, 'printing a, b, c e s', a, b, c, s 
      !print*, 'printing tan', tan(0.5*s), tan(0.5*(s-a)), tan(0.5*(s-b)),
      !tan(0.5*(s-c))

      tanqe = sqrt(tan(0.5*s)*tan(0.5*(s-a))*tan(0.5*(s-b))*tan(0.5*(s-c)))
      sphericalTriangleArea = 4.*atan(tanqe)*radius**2

      !print*, 'printing a, b, c, s e tanqe', a, b, c, s, tanqe

   end function sphericalTriangleArea

   real function sphere_distance(lat1, lon1, lat2, lon2, radius)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on
   ! a sphere with given radius.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius

      real (kind=RKIND) :: arg1

      arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
                   cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
      sphere_distance = 2.*radius*asin(arg1)

   end function sphere_distance


end module ocn_ke_triskCV
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
