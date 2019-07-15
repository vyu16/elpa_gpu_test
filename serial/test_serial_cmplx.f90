! This program tests ELPA.
program test_serial_cmplx

   use ELPA

   implicit none

   integer, parameter :: dp = selected_real_kind(15,300)

   character(len=10) :: arg1
   character(len=10) :: arg2
   character(len=10) :: arg3
   character(len=10) :: arg4
   character(len=10) :: arg5

   integer :: ierr
   integer :: blk
   integer :: n_basis
   integer :: n_states
   integer :: method
   integer :: cpu
   integer :: n
   integer :: i

   real(dp) :: err1
   real(dp) :: err2
   complex(dp) :: aux

   integer, allocatable :: seed(:)

   complex(dp), allocatable :: mat(:,:)
   complex(dp), allocatable :: tmp(:,:)
   real(dp), allocatable :: help(:,:)
   complex(dp), allocatable :: evec(:,:)
   real(dp), allocatable :: eval(:)

   class(elpa_t), pointer :: eh

   complex(dp), external :: zdotc

   complex(dp), parameter :: zero = (0.0_dp,0.0_dp)
   complex(dp), parameter :: one = (1.0_dp,0.0_dp)

   ! Read command line arguments
   if(COMMAND_ARGUMENT_COUNT() == 5) then
      call GET_COMMAND_ARGUMENT(1,arg1)
      call GET_COMMAND_ARGUMENT(2,arg2)
      call GET_COMMAND_ARGUMENT(3,arg3)
      call GET_COMMAND_ARGUMENT(4,arg4)
      call GET_COMMAND_ARGUMENT(5,arg5)

      read(arg1,*) n_basis
      if(n_basis <= 0) then
         n_basis = 1000
      end if

      read(arg2,*) n_states
      if(n_states < 0 .or. n_states > n_basis) then
         n_states = n_basis
      end if

      read(arg3,*) blk
      read(arg4,*) method
      read(arg5,*) cpu
   else
      write(*,"(2X,A)") "################################################"
      write(*,"(2X,A)") "##  Wrong number of command line arguments!!  ##"
      write(*,"(2X,A)") "##  Arg#1: Size of test matrix.               ##"
      write(*,"(2X,A)") "##  Arg#2: Number of eigenvectors to compute. ##"
      write(*,"(2X,A)") "##  Arg#3: Block size                         ##"
      write(*,"(2X,A)") "##  Arg#4: 1 = ELPA 1-stage                   ##"
      write(*,"(2X,A)") "##         2 = ELPA 2-stage                   ##"
      write(*,"(2X,A)") "##  Arg#5: 1 = CPU                            ##"
      write(*,"(2X,A)") "##         2 = GPU                            ##"
      write(*,"(2X,A)") "################################################"
      stop
   end if

   ! Generate a random matrix
   call random_seed(size=n)

   allocate(seed(n))
   allocate(mat(n_basis,n_basis))
   allocate(tmp(n_basis,n_basis))
   allocate(help(n_basis,n_basis))

   seed = 12345678

   call random_seed(put=seed)
   call random_number(help)

   ! Symmetrize test matrix
   tmp = help+(0.0_dp,1.0_dp)*help
   mat = tmp+transpose(conjg(tmp))
   tmp = mat

   deallocate(help)
   allocate(evec(n_basis,n_basis))
   allocate(eval(n_basis))

   write(*,*)
   write(*,"(2X,A)") "Test matrix generated"
   write(*,*)

   ! Initialize ELPA
   ierr = elpa_init(20180525)
   eh => elpa_allocate()

   call eh%set("na",n_basis,ierr)
   call eh%set("nev",n_states,ierr)
   call eh%set("nblk",blk,ierr)
   call eh%set("local_nrows",n_basis,ierr)
   call eh%set("local_ncols",n_basis,ierr)
   call eh%set("timings",1,ierr)

   ierr = eh%setup()

   if(ierr /= 0) then
      write(*,"(2X,A)") "Error: setup"
      stop
   end if

   if(method == 1) then
      call eh%set("solver",ELPA_SOLVER_1STAGE,ierr)
   else
      call eh%set("solver",ELPA_SOLVER_2STAGE,ierr)
   end if

   if(cpu == 1) then
      call eh%set("gpu",0,ierr)
      call eh%set("complex_kernel",ELPA_2STAGE_COMPLEX_GENERIC,ierr)
   else
      call eh%set("gpu",1,ierr)
      call eh%set("complex_kernel",ELPA_2STAGE_COMPLEX_GPU,ierr)
   end if

   ! Solve
   call eh%eigenvectors(mat,eval,evec,ierr)

   if(ierr /= 0) then
      write(*,"(2X,A)") "Error: solve"
      stop
   end if

   write(*,"(2X,A)") "ELPA solver finished"

!   call eh%print_times()

   ! Finalize ELPA
   call elpa_deallocate(eh)

   nullify(eh)

   ! Check A C - lambda C
   call zgemm("N","N",n_basis,n_states,n_basis,one,tmp,n_basis,evec,n_basis,&
        zero,mat,n_basis)

   tmp = evec

   do i = 1,n_states
      tmp(:,i) = tmp(:,i)*eval(i)
   end do

   tmp = tmp-mat
   err1 = 0.0_dp

   do i = 1,n_states
      aux = zdotc(n_basis,tmp(:,i),1,tmp(:,i),1)

      err1 = max(err1,sqrt(real(aux,kind=dp)))
   enddo

   ! Check I - C^T C
   tmp = zero

   call zlaset("U",n_states,n_states,zero,one,tmp,n_basis)
   call zherk("U","C",n_states,n_basis,-one,evec,n_basis,one,tmp,n_basis)

   err2 = maxval(abs(tmp))

   write(*,"(2X,A,E10.2,A,E10.2)") "| Error :",err1,";",err2
   write(*,*)

   deallocate(mat)
   deallocate(eval)
   deallocate(evec)
   deallocate(tmp)

end program
