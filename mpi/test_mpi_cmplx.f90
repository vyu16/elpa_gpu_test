! This program tests ELPA.
program test_mpi_cmplx

   use ELPA
   use MPI

   implicit none

   integer, parameter :: dp = selected_real_kind(15,300)

   character(len=10) :: arg1
   character(len=10) :: arg2
   character(len=10) :: arg3
   character(len=10) :: arg4
   character(len=10) :: arg5

   integer :: n_proc
   integer :: nprow
   integer :: npcol
   integer :: myid
   integer :: myprow
   integer :: mypcol
   integer :: comm
   integer :: ierr
   integer :: blk
   integer :: ctxt
   integer :: desc(9)
   integer :: n_basis
   integer :: n_states
   integer :: method
   integer :: cpu
   integer :: nlrow
   integer :: nlcol
   integer :: ldm
   integer :: n
   integer :: i

   real(dp) :: t1
   real(dp) :: t2
   real(dp) :: myerr
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

   integer, external :: numroc

   complex(dp), parameter :: zero = (0.0_dp,0.0_dp)
   complex(dp), parameter :: one = (1.0_dp,0.0_dp)

   ! Initialize MPI
   call MPI_Init(ierr)
   comm = MPI_COMM_WORLD
   call MPI_Comm_size(comm,n_proc,ierr)
   call MPI_Comm_rank(comm,myid,ierr)

   ! Read command line arguments
   if(command_argument_count() == 5) then
      call get_command_argument(1,arg1)
      call get_command_argument(2,arg2)
      call get_command_argument(3,arg3)
      call get_command_argument(4,arg4)
      call get_command_argument(5,arg5)

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
      if(myid == 0) then
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
         call MPI_Abort(comm,0,ierr)
         stop
      end if
   end if

   ! Set up square-like processor grid
   do npcol = nint(sqrt(real(n_proc))),2,-1
      if(mod(n_proc,npcol) == 0) exit
   end do
   nprow = n_proc/npcol

   ! Set up BLACS
   ctxt = comm

   call BLACS_Gridinit(ctxt,"r",nprow,npcol)
   call BLACS_Gridinfo(ctxt,nprow,npcol,myprow,mypcol)

   nlrow = numroc(n_basis,blk,myprow,0,nprow)
   nlcol = numroc(n_basis,blk,mypcol,0,npcol)

   ldm = max(nlrow,1)

   call descinit(desc,n_basis,n_basis,blk,blk,0,0,ctxt,ldm,ierr)

   ! Generate a random matrix
   call random_seed(size=n)

   allocate(seed(n))
   allocate(mat(nlrow,nlcol))
   allocate(tmp(nlrow,nlcol))
   allocate(help(nlrow,nlcol))

   seed(:) = myid+1

   call random_seed(put=seed)
   call random_number(help)

   ! Symmetrize test matrix
   mat(:,:) = help+(0.0_dp,1.0_dp)*help
   tmp(:,:) = mat

   call pztranc(n_basis,n_basis,one,tmp,1,1,desc,one,mat,1,1,desc)

   ! Save test matrix
   tmp(:,:) = mat

   deallocate(help)
   allocate(evec(nlrow,nlcol))
   allocate(eval(n_basis))

   if(myid == 0) then
      write(*,*)
      write(*,"(2X,A)") "Complex test matrix generated"
      write(*,"(2X,A,I10)") "| Matrix size  :",n_basis
      write(*,"(2X,A,I10)") "| Eigenvectors :",n_states
      write(*,"(2X,A,I10)") "| Block size   :",blk
      write(*,"(2X,A,I10)") "| MPI tasks    :",n_proc
      write(*,*)
   end if

   ! Initialize ELPA
   ierr = elpa_init(20180525)
   eh => elpa_allocate(ierr)

   call eh%set("na",n_basis,ierr)
   call eh%set("nev",n_states,ierr)
   call eh%set("nblk",blk,ierr)
   call eh%set("local_nrows",nlrow,ierr)
   call eh%set("local_ncols",nlcol,ierr)
   call eh%set("mpi_comm_parent",comm,ierr)
   call eh%set("process_row",myprow,ierr)
   call eh%set("process_col",mypcol,ierr)
   call eh%set("timings",1,ierr)

   ierr = eh%setup()

   if(ierr /= 0) then
      write(*,"(2X,A)") "Error: setup"
      call MPI_Abort(comm,0,ierr)
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

   t1 = MPI_Wtime()

   ! Solve
   call eh%eigenvectors(mat,eval,evec,ierr)

   t2 = MPI_Wtime()

   if(ierr /= 0) then
      write(*,"(2X,A)") "Error: solve"
      call MPI_Abort(comm,0,ierr)
      stop
   end if

   if(myid == 0) then
      write(*,"(2X,A)") "ELPA solver finished"
      write(*,"(2X,A,F10.3,A)") "| Time  :",t2-t1,"s"

      call eh%print_times()
   end if

   ! Finalize ELPA
   call elpa_deallocate(eh,ierr)

   nullify(eh)

   if(.false.) then
      ! Check A C - lambda C
      call pzgemm("N","N",n_basis,n_states,n_basis,one,tmp,1,1,desc,evec,1,1,desc,&
           zero,mat,1,1,desc)

      tmp = evec

      do i = 1,n_states
         aux = eval(i)

         call pzscal(n_basis,aux,tmp,1,i,desc,1)
      end do

      tmp = tmp-mat
      myerr = 0.0_dp

      do i = 1,n_states
         call pzdotc(n_basis,aux,tmp,1,i,desc,1,tmp,1,i,desc,1)

         myerr = max(myerr,sqrt(real(aux,kind=dp)))
      end do

      call MPI_Reduce(myerr,err1,1,MPI_REAL8,MPI_MAX,0,comm,ierr)

      ! Check I - C^T C
      tmp = zero

      call pzlaset("U",n_states,n_states,zero,one,tmp,1,1,desc)
      call pzherk("U","C",n_states,n_basis,-one,evec,1,1,desc,one,tmp,1,1,desc)

      myerr = maxval(abs(tmp))

      call MPI_Reduce(myerr,err2,1,MPI_REAL8,MPI_MAX,0,comm,ierr)

      if(myid == 0) then
         write(*,"(2X,A,E10.2,A,E10.2)") "| Error :",err1,";",err2
         if(err1 > 1.e-9_dp .or. err2 > 1.e-11_dp) then
            write(*,"(2X,A)") "Failed!!"
         end if
         write(*,*)
      end if
   end if

   deallocate(mat)
   deallocate(eval)
   deallocate(evec)
   deallocate(tmp)
   deallocate(seed)

   call BLACS_Gridexit(ctxt)
   call BLACS_Exit(1)
   call MPI_Finalize(ierr)

end program
