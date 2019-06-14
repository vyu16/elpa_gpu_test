! This program tests ELPA.
program test_elpa

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
   integer :: info
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

   integer, allocatable :: seed(:)

   real(dp), allocatable :: mat(:,:)
   real(dp), allocatable :: tmp(:,:)
   real(dp), allocatable :: evec(:,:)
   real(dp), allocatable :: eval(:)

   class(elpa_t), pointer :: eh

   integer, external :: numroc

   ! Initialize MPI
   call MPI_Init(ierr)
   comm = MPI_COMM_WORLD
   call MPI_Comm_size(comm,n_proc,ierr)
   call MPI_Comm_rank(comm,myid,ierr)

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

   call descinit(desc,n_basis,n_basis,blk,blk,0,0,ctxt,ldm,info)

   ! Generate a random matrix
   call random_seed(size=n)

   allocate(seed(n))
   allocate(mat(nlrow,nlcol))
   allocate(tmp(nlrow,nlcol))

   seed = myid

   call random_seed(put=seed)
   call random_number(mat)

   ! Symmetrize test matrix
   tmp = mat

   call pdtran(n_basis,n_basis,1.0_dp,tmp,1,1,desc,1.0_dp,mat,1,1,desc)

   ! Save test matrix
   tmp = mat

   allocate(evec(nlrow,nlcol))
   allocate(eval(n_basis))

   if(myid == 0) then
      write(*,*)
      write(*,"(2X,A)") "Test matrix generated"
      write(*,*)
   end if

   ! Initialize ELPA
   ierr = elpa_init(20180525)
   eh => elpa_allocate()

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
   end if

   if(method == 1) then
      call eh%set("solver",ELPA_SOLVER_1STAGE,ierr)
   else
      call eh%set("solver",ELPA_SOLVER_2STAGE,ierr)
   end if

   if(cpu == 1) then
      call eh%set("gpu",0,ierr)
      call eh%set("real_kernel",ELPA_2STAGE_REAL_GENERIC,ierr)
      call eh%set("complex_kernel",ELPA_2STAGE_COMPLEX_GENERIC,ierr)
   else
      call eh%set("gpu",1,ierr)
      call eh%set("real_kernel",ELPA_2STAGE_REAL_GPU,ierr)
      call eh%set("complex_kernel",ELPA_2STAGE_COMPLEX_GPU,ierr)
   end if

   t1 = MPI_Wtime()

   ! Solve
   call eh%eigenvectors(mat,eval,evec,ierr)

   t2 = MPI_Wtime()

   if(ierr /= 0) then
      write(*,"(2X,A)") "Error: solve"
   end if

   if(myid == 0) then
      write(*,"(2X,A)") "ELPA solver finished"
      write(*,"(2X,A,F10.3,A)") "| Time  :",t2-t1,"s"

!      call eh%print_times()
   end if

   ! Finalize ELPA
   call elpa_deallocate(eh)

   nullify(eh)

   ! Check A C - lambda C
   call pdgemm("N","N",n_basis,n_states,n_basis,1.0_dp,tmp,1,1,desc,evec,1,1,&
        desc,0.0_dp,mat,1,1,desc)

   tmp = evec

   do i = 1,n_states
      call pdscal(n_basis,eval(i),tmp,1,i,desc,1)
   end do

   tmp = tmp-mat
   myerr = 0.0_dp

   do i = 1,n_states
      call pdnrm2(n_basis,err1,tmp,1,i,desc,1)

      myerr = max(myerr,err1)
   enddo

   call MPI_Reduce(myerr,err1,1,MPI_REAL8,MPI_MAX,0,comm,ierr)

   ! Check I - C^T C
   tmp = 0.0_dp

   call pdlaset("U",n_states,n_states,0.0_dp,1.0_dp,tmp,1,1,desc)
   call pdsyrk("U","T",n_states,n_basis,-1.0_dp,evec,1,1,desc,1.0_dp,tmp,1,1,&
        desc)

   myerr = maxval(abs(tmp))

   call MPI_Reduce(myerr,err2,1,MPI_REAL8,MPI_MAX,0,comm,ierr)

   if(myid == 0) then
      write(*,"(2X,A,E10.2,A,E10.2)") "| Error :",err1,";",err2
      write(*,*)
   end if

   deallocate(mat)
   deallocate(eval)
   deallocate(evec)
   deallocate(tmp)

   call BLACS_Gridexit(ctxt)
   call BLACS_Exit(1)
   call MPI_Finalize(ierr)

end program