!----------------------------------------------------------------------------------------------------------------------------------------   
! OpenCFD Ver 1.9,  3-D compressible Navier-Stokes Finite difference Solver 
! Copyright by LI Xinliang, LHD, Institute of Mechanics, CAS, Email: lixl@imech.ac.cn
!  
! The default code is double precision computation
! If you want to use SINGLE PRECISION computation, you can change   "OCFD_REAL_KIND=8"  to "OCFD_REAL_KIND=4" ,
! and  "OCFD_DATA_TYPE=MPI_DOUBLE_PRECISION" to "OCFD_DATA_TYPE=MPI_REAL" in the file OpenCFD.h 
!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------------- 
   program main
      Use OCFD_flow_parameters
      implicit none

      integer,parameter:: IBUFFER_SIZE=1000000    
      real(kind=8)  BUFFER_MPI(IBUFFER_SIZE)    ! Buffer for MPI  message transfer (used by MPI_Bsend)
      integer   ierr, status(MPI_status_size)

!------------------------------------------------
       call mpi_init(ierr)
       call mpi_comm_rank(MPI_COMM_WORLD,my_id,ierr)
       call mpi_comm_size(MPI_COMM_WORLD,np_size,ierr)
       call MPI_BUFFER_ATTACH(BUFFER_MPI,8*IBUFFER_SIZE,ierr)

!c npx0,npy0,npz0 are parallel partations in x,y and z directions    
        call read_parameters
        call part
        call weno
        call mpi_finalize(ierr)

   end
!------------------------------------------------------------------------------------------------
