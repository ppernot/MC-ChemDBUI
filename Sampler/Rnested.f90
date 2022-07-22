! Generation of Monte Carlo samples from Dirichlet-type distributions
! as defined in Plessis, Carrasco & Pernot 2010
!
! This file contains two main modules:
! 1/ Dirichlet: set of basic (i.e. non-nested) distributions and utilities
! 2/ eagle: gestion of data trees derived and adapted from the code of Metcalf
! The main program reads the number of samples and the expression of the distribution and
! outputs an ascii file with the sample points.
!
! For random numbers generation, the code should be linked with ranlib (gengam)
! e.g. ifort -O3 nested.f90 /home/lib/ranlib.a -o nested.x

! 2014-01-03 (PP): added constriants to 3*sigma for unbounded distribs (Dirg, Logn)

module Dirichlet
! Set of basic (i.e. non-nested) distributions and utilities
contains

 function Mlogn(alpha,beta) result (rand_tab)
  ! Multivatiate lognormal with independent variables
  implicit none
  real, dimension(:),intent(in) :: alpha, beta
  real, dimension(size(alpha))  :: rand_tab
  real                          :: rn, snorm
  logical                       :: ok
  integer                       :: i

  do i= 1, size(alpha)
    ok = .FALSE.
    do while (.NOT. ok)
       rn = snorm()
       if(rn .GT. -3 .AND. rn .LT. 3) ok=.TRUE.
    enddo
    rand_tab(i) = exp(log(alpha(i)) + log(beta(i))*rn)
  enddo

 end function Mlogn

 function Logn(alpha) result (rand_tab)

  implicit none
  real, dimension(2),intent(in)   :: alpha
  real                            :: rand_tab
  real                            :: rn, snorm
  logical                         :: ok

  ok = .FALSE.
  do while (.NOT. ok)
     rn = snorm()
     if(rn .GT. -3 .AND. rn .LT. 3) ok=.TRUE.
  enddo
  rand_tab=exp(log(alpha(1)) + log(alpha(2))*rn)

  end function Logn

  function Logu(alpha) result (rand_tab)

  implicit none
  real, dimension(2),intent(in)   :: alpha
  real                            :: rand_tab
  real                            :: genunf

  rand_tab=exp(genunf(log(alpha(1)),log(alpha(2))))

  end function Logu

  function Norm(alpha) result (rand_tab)

  implicit none
  real, dimension(2),intent(in)   :: alpha
  real                            :: rand_tab
  real                            :: rn, snorm
  logical                         :: ok

  ok = .FALSE.
  do while (.NOT. ok)
     rn = snorm()
     if(rn .GT. -3 .AND. rn .LT. 3) ok=.TRUE.
  enddo
  rand_tab=alpha(1)+alpha(2)*rn

  end function Norm

  function Nort(alpha,beta) result (rand_tab)

  implicit none
  real, dimension(2),intent(in)   :: alpha, beta
  real                            :: rand_tab
  real                            :: snorm

  do
    rand_tab=alpha(1)+alpha(2)*snorm()
    if ((rand_tab.ge.min(beta(1),beta(2))).and.(rand_tab.le.max(beta(1),beta(2)))) exit
  enddo

  end function Nort

  function Unif(alpha) result (rand_tab)

  implicit none
  real, dimension(2),intent(in)   :: alpha
  real                            :: rand_tab
  real                            :: genunf

  rand_tab=genunf(alpha(1),alpha(2))

  end function Unif

  function Drat(alpha,beta) result (rand_tab)
  ! Generate two variables with unit sum and variable ratio
  ! x1/x2= alpha */: beta (beta is uncertainty factor >1)
    implicit none
    real, intent(in)              :: alpha, beta
    real, dimension(2)            :: rand_tab
    real                          :: u, snorm

    ! Generate random ratio with positivity constraint
    u   = exp(log(alpha) + log(beta)*snorm())
    rand_tab(2)=1/(1+u)
    rand_tab(1)=1-rand_tab(2)

  end function Drat

  function Disc(alpha) result (rand_tab)
  ! Discrete sampling according to probas alpha/sum(alpha)
    implicit none
    real, dimension(:),intent(in) :: alpha
    real, dimension(size(alpha))  :: rand_tab
    real                          :: u, sa, vmin, vmax,genunf
    integer                       :: i

    rand_tab=0
    u   = genunf(0.0,1.0)
    sa  =sum(alpha)
    vmax=0
    do i= 1, size(alpha)
      vmin=vmax
      vmax=vmax+alpha(i)/sa
      if((u-vmin)*(u-vmax).le.0) then
         rand_tab(i)= 1
         exit
      endif
    enddo

  end function Disc

  function Dirw(alpha,beta) result (rand_tab)
  ! Generalized Dirichlet distribution according to Wong1998
  ! tolerating invariables in the distribution
    implicit none
    real, dimension(:),intent(in) :: alpha, beta
    real, dimension(size(alpha))  :: rand_tab
    real                          :: genbet, a, b
    integer                       :: i

    a=alpha(1)
    b=beta(1)
    rand_tab(1)=genbet(a,b)
    do i= 2, size(alpha)-1
       a=alpha(i)
       b=beta(i)
       rand_tab(i)=genbet(a,b)*(1.0e0-sum(rand_tab(1:i-1))) ! Generate Beta deviates
    enddo
    rand_tab(size(alpha))=1.0e0-sum(rand_tab(1:i-1))

  end function Dirw

  function Dirz(alpha,beta) result (rand_tab)
  ! Generalized Dirichlet distribution according to Wong1998
  ! tolerating invariables in the distribution
  ! Differs from Dirw: the parameters are the mean and stdv
    implicit none
    real, dimension(:),intent(in) :: alpha, beta
    real, dimension(size(alpha))  :: rand_tab, rand_tabc
    real                          :: genbet, a, b, si, p1, p2, pm, pv
    integer                       :: i, ic
    real, parameter               :: beta_min=1.0e-3, & ! lower limit of significant uncertainty
                                     sf=1.0e-0          ! scaling factor

    rand_tab=alpha
    do
      rand_tabc=0e0
      ic=0
      si=0e0
      pm=1.0e0
      pv=1.0e0
      do i= 1, size(alpha)
        if (beta(i).LE.beta_min) then
          if ((i.EQ.size(alpha)).AND.(beta(i).LT.0e0)) then
            ic=ic+1
            rand_tabc(ic)=sf*(1.0e0-si)-sum(rand_tabc)
          else
            si=si+alpha(i)
          endif
        else
          ic=ic+1
          p1=alpha(i)*sf/pm
          p2=(beta(i)**2+alpha(i)**2)*sf**2/pv
          a=p1*(p1-p2)/(p2-p1**2)
          b=(1.0e0-p1)*(p1-p2)/(p2-p1**2)
          pm=pm*b/(a+b)
          pv=pv*b/(a+b)*(b+1.0e0)/(a+b+1.0e0)
          rand_tabc(ic)=genbet(a,b)*(1.0e0-sum(rand_tabc(1:ic-1))) ! Generate Beta deviates
        endif
      enddo
      if (rand_tabc(ic).GE.0e0) exit
    enddo
    if (ic.NE.0) then
      if (si.GE.1.0e0) stop 'Nominal values too large'
      ic=0
      rand_tabc=rand_tabc/sum(rand_tabc)*(1.0e0-si)
      do i= 1, size(alpha)
        if ((beta(i).LE.beta_min).AND.((beta(i).GE.0e0).OR.(i.NE.size(alpha)))) cycle
        ic=ic+1
        rand_tab(i)=rand_tabc(ic)
      enddo
    endif

  end function Dirz

  function Dirg(alpha,beta) result (rand_tab)
  ! Generalized Dirichlet distribution according to Lingwan REF???
    implicit none
    real, dimension(:),intent(in) :: alpha, beta
    real, dimension(size(alpha))  :: rand_tab
    real                          :: gengam, a, b
    integer                       :: i, count
    logical                       :: ok

    ok = .FALSE.
    count = 0
    do while (.NOT. ok)

       do i= 1, size(alpha)
          a= alpha(i)/beta(i)**2
          b=(alpha(i)/beta(i))**2
          rand_tab(i)= gengam(a,b) ! Generate Gamma deviates
       enddo
       rand_tab=rand_tab/sum(rand_tab)

       ! Control of stray runs by marginals
       ok = .TRUE.
       do i= 1, size(alpha)
          if(rand_tab(i) .LT. alpha(i)-3*beta(i) .OR. &
             rand_tab(i) .GT. alpha(i)+3*beta(i) ) then
             ok = .FALSE.
             count = count +1
             if (count .GT.50) stop 'Problem in Dirg'
             exit
          endif
       enddo

    enddo

  end function Dirg

  function Diri(alpha,beta) result (rand_tab)
  ! Dirichlet distribution (cf. Gelman)
    implicit none
    real, dimension(:),intent(in) :: alpha, beta
    real, dimension(size(alpha))  :: rand_tab
    real                          :: sgamma, a
    integer                       :: i

    do i=  1, size(alpha)
      a=alpha(i)*beta(1)
      rand_tab(i)= sgamma(a) ! Generate Gamma deviates
    enddo
    rand_tab=rand_tab/sum(rand_tab)

  end function Diri

  function Diun(n) result (rand_tab)
  ! Uniform Dirichlet
    implicit none
    integer, intent(in)  :: n
    real, dimension(n) :: rand_tab, ones

    ones=1
    rand_tab= Dirg(ones,ones)

  end function Diun

  function Dior(alpha) result (rand_tab)
  ! Ordered Uniform Dirichlet
    implicit none
    real, dimension(:)             :: alpha
    integer                        :: n, i
    real, dimension(size(alpha))   :: rand_tab, wrk

    n=size(alpha)
    wrk=Diun(n)
    wrk=sort_down(wrk)
    do i= 1, n
       rand_tab(i)=wrk(nint(alpha(i)))
    enddo

  end function Dior

  function Diut(alpha,beta) result (rand_tab)
  ! Intervals-constrained Uniform Dirichlet
  ! Algorithm of Fang & Yang, Stat. Prob. Letters 46:113-120 (2000)
    implicit none
    real, dimension(:)          :: alpha, beta
    real, dimension(size(alpha)):: ai, bi
    real, dimension(size(alpha)):: rand_tab
    real                        :: a, b, u, dk, delk, phik
    real                        :: genunf, b1, b2
    integer                     :: n, i, k

    ai=alpha
    bi=beta
    n=size(alpha)

    ! O/ Consistency check:
    !    Sum of lower bounds cannot be larger than one
    !    Sum of upper bounds cannot be smaller than one
    a=sum(ai)
    b=sum(bi)
    if (a.GE.1 .OR. b.LE.1) stop 'Diut: bounds consistency problem'

    ! 1/ Remove spurious constraints
    do i= 1, n
       b1=max(ai(i),bi(i)+1-b)
       b2=min(bi(i),ai(i)+1-a)
       ai(i)=b1
       bi(i)=b2
    enddo
    a=sum(ai)
    b=sum(bi)

    ! 2/ Generate variables by uniform sampling
    !    with successive bounds constraints
    delk=1
    do k= n, 2, -1
       u   = genunf(0.0,1.0)
       dk  = max(ai(k)/delk,1.0-sum(bi(1:k-1))/delk)
       phik= min(bi(k)/delk,1.0-sum(ai(1:k-1))/delk)
       rand_tab(k)=G_Fang(u,dk,phik,delk,k-1)
       delk= 1-sum(rand_tab(k:n))
    enddo
    rand_tab(1)=delk

    contains
      real function G_Fang(u,d,b,c,k)
        ! G function of Fang & Yang, Stat. Prob. Letters 46:113-120 (2000)
        implicit none
        real    :: u,d,b,c
        integer :: k

        G_Fang= c * ( 1 - ( u*(1-b)**k + (1-u)*(1-d)**k )**(1.d0/k) )

      end function G_Fang

  end function Diut


  function simplex_3D(x) result(y)
  ! Convert x,y,z into projection onto 3-simplex
    implicit none
    real, dimension(3) :: x
    real, dimension(2) :: y
    real               :: rac2, rac6

    rac2=sqrt(2.0)
    rac6=sqrt(6.0)

    y(1)=-(x(1)-1./3.)/rac2+(x(2)-1./3.)/rac2
    y(2)=-(x(1)-1./3.)/rac6-(x(2)-1./3.)/rac6+2*(x(3)-1./3.)/rac6

  end function simplex_3D

  function sort_down(x) result(y)
  ! Sort table by decreasing order
    implicit none
    real, dimension(:)       :: x
    real, dimension(size(x)) :: y
    real :: swap
    integer :: change, i, n

    n=size(x)

    change=1
    do while (change.eq.1)
      change=0
      do i= 1, n-1
        if(x(i) .lt. x(i+1)) then
          change=1
          swap=x(i+1)
          x(i+1)=x(i)
          x(i)=swap
        endif
      enddo
    enddo
    y=x

  end function sort_down

end module Dirichlet


! Simplified version of Metcalf's code available at
! http://www.nag.co.uk/nagware/Examples/trees.f90

! (c) Copyright Michael Metcalf and John Reid, 1992. This file may be
! freely used and copied for educational purposes provided this notice
! remains attached. Extracted from "Fortran 90 Explained" Oxford
! University Press (Oxford and New York), ISBN 0-19-853772-7.
!
module ddl
! Contains a type definition that is to become a component of the node data type.
   type user_type
      character(8)     :: layout
      integer, pointer :: i_field(:)
      real, pointer    :: f_field(:)
   end type user_type

end module ddl

module eagle    ! Soars above the rest
   use ddl

! Stong typing imposed
   implicit none
!
! Only the data type, the subroutine interfaces, the length of the character
! component, the user type are public
   private
   public finish, new_node, start, user_type, max_char, &            ! Metcalf's routines
          rng_tree, num_leaves, build_tree, max_leaves, chain_size,& ! My routines
          nestdiri                                                   ! My routines

! Global constants
   integer, parameter     :: max_char  = 32  ! length of character component

! My globals
   integer, parameter     :: max_leaves= 50  ! max number of children for a node
   integer, parameter     :: chain_size= 2048 ! max size of char chain to process
   integer, parameter     :: max_nodes = 500 ! max number of nodes in tree

! Define a pointer data type
   type ptr
      type(data), pointer :: pp
   end type ptr
!
! Define the data structure holding the state of a tree
   type state
      type(data), pointer :: current, parent
      integer             :: count, last_code
   end type state
!
! Define the basic node type
   type, public :: data
   private
      integer                     index, amount(4), running_index
      character(max_char)         header
      integer, pointer         :: j(:), link(:)
      real, pointer            :: y(:)
      type(ptr), pointer       :: p(:)
      type(user_type), pointer :: user(:)
      type(data), pointer      :: back, reference, ref_back
      type(state), pointer     :: own_state
   end type data
!
! Some global module variables
   type(data), pointer    :: current, parent
   integer                :: count, last_code
   type(state), pointer   :: tree_state        ! tree state variable
!
! The module procedures
contains

   subroutine add_user(node_out, node_in)
!
! To make the assigment of the user data to a node. Knowlege of
! the actual structure of user_type is required.
      type(user_type), intent(out) :: node_out(:)
      type(user_type), intent(in)  :: node_in(:)
      integer loop
!
      do loop = 1, size(node_in)
         node_out(loop)%layout = node_in(loop)%layout
         allocate(node_out(loop)%i_field(size(node_in(loop)%i_field)))
         allocate(node_out(loop)%f_field(size(node_in(loop)%f_field)))
         node_out(loop)%i_field = node_in(loop)%i_field
         node_out(loop)%f_field = node_in(loop)%f_field
      end do
   end subroutine add_user

   subroutine check_state(tree)
!
! To check that this is the same tree as on the last call, and to
! switch the state variable if not.
      type(data), intent(in) :: tree
!
      if (.not.associated(tree%own_state, tree_state)) then
         tree_state => tree%own_state
         count =  tree_state%count
         last_code = tree_state%last_code
         current => tree_state%current
         parent => tree_state%parent
      end if
   end subroutine check_state

   recursive subroutine finish (tree)
!
! Traverse a complete tree or subtree, deallocating all storage
! (except the state variable).
      type(data), pointer :: tree
      integer loop
!
      do loop = 1, size(tree%p)
         call finish (tree%p(loop)%pp)
      end do
      do loop = 1, size(tree%user)
         deallocate(tree%user(loop)%i_field, tree%user(loop)%f_field)
      end do
      deallocate(tree%j, tree%y, tree%user, tree%p, tree%link)
      if (associated(tree%ref_back)) nullify(tree%ref_back%reference)
      deallocate(tree)
   end subroutine finish

   recursive subroutine new_data(node, number, node_link, name, real_data,     &
                                 integer_data, user_data, links)
!
! Add a new node to the tree
      type(data), pointer                :: node
      integer, intent(in)                :: number
      character(*), optional, intent(in) :: name
      real, optional, intent(in)         :: real_data(:)
      type(user_type), optional, intent(in) :: user_data(:)
      integer, optional, intent(in)      :: integer_data(:), links(:)
      type(data), pointer    :: other, first
      integer                :: loop, control, node_link
!
! Save the old counter
      control = count
!
! Test whether a node is already there, if it is go down one level
      if (associated(node)) then
         do loop = 1, size(node%p)
            call new_data(node%p(loop)%pp, number, node%link(loop), name,      &
                         real_data = real_data,                                &
                         integer_data = integer_data, user_data = user_data,   &
                         links = links)
!
! Modify back pointer if depends on higher layer
            if (loop == 1 .and. size(node%p) > 1 ) then
               first => node%p(1)%pp
               if(associated(first)) first => first%back    ! Precaution
            end if
!
! Have we added the node? If so reset the back pointer if necessary and exit.
            if (count > control) then
               if (loop > 1 ) then
                  other => node%p(loop)%pp
                  other%back => first
               end if
               exit
            endif
         end do
!
! If it isn't there, add it
      elseif (count == 0 .or. number == node_link) then
         count = count + 1
         allocate(node)
         node%own_state => tree_state
!
! Set the back pointer to a preliminary value
         if (count == 1) then
            nullify(node%back)
         else
            node%back => current
         end if
         node%index = number
         node%running_index = count
         nullify(node%reference, node%ref_back)
!
! Add the data where present
         if (present(name)) then
            node%header = name
         else
            node%header = ' '
         end if
         if (present(integer_data)) then
            allocate(node%j(size(integer_data)))
            node%j = integer_data
            node%amount(1) = size(integer_data)
         else
            allocate(node%j(0))
            node%amount(1) = 0
         end if
         if (present(real_data)) then
            allocate(node%y(size(real_data)))
            node%y = real_data
            node%amount(2) = size(real_data)
         else
            allocate(node%y(0))
            node%amount(2) = 0
         end if
         if (present(user_data)) then
            allocate(node%user(size(user_data)))
            call add_user(node%user, user_data)
            node%amount(4) = size(user_data)
         else
            allocate(node%user(0))
            node%amount(4) = 0
         end if
!
! Add links if present and valid
         if (present(links)) then
            if (size(links) > 0) then
               if (all(links > count)) then
                  allocate(node%p(size(links)), node%link(size(links)))
                  current => node
                  node%link = links
                  node%amount(3) = size(links)
                  do loop = 1, size(links)
                     nullify(node%p(loop)%pp)
                  end do
               else
                  print *,  'attempt to link node',  count,   'to higher-level &
                          &(node one of:', links, ')'
                  stop 'link error'
               end if
            else
               allocate(node%p(0), node%link(0))
               node%amount(3) = 0
            end if
         else
            allocate(node%p(0), node%link(0))
            node%amount(3) = 0
         end if
      end if
   end subroutine new_data

   subroutine new_node(tree, number, name, real_data, integer_data, user_data, &
                       links)
!
! To add a node efficiently and to check that a node actually gets added
      type(data), pointer                :: tree
      integer, intent(in)                :: number
      character(*), optional, intent(in) :: name
      real, optional, intent(in)         :: real_data(:)
      type(user_type), optional, intent(in) :: user_data(:)
      integer, optional, intent(in)      :: integer_data(:), links(:)
      type(data), pointer    :: bottom
      integer                :: control, dummy_link = 0
!
      if (count /= 0) call check_state(tree)
!
! First try to attach to last node added,
      control = count
      if (associated(current)) then
         if (size(current%p) /= 0) then
            bottom => current
            call new_data(bottom, number, dummy_link, name,                    &
                                        real_data = real_data,                 &
                   integer_data = integer_data, user_data = user_data,         &
                   links = links)
         end if
      endif
      if (count == control) then
!
! otherwise, do complete search.
         call new_data(tree, number, dummy_link, name,                         &
                                   real_data = real_data,                      &
                   integer_data = integer_data, user_data = user_data,         &
                   links = links)
         if (count == control) print *, 'new node', number,                    &
                                  ' could not be added to closed tree'
      endif
      call update_state(tree)
   end subroutine new_node

   subroutine start(tree)
!
! THE CALL TO START FOR A GIVEN TREE MUST IMMEDIATELY PRECEDE
! THE FIRST CALL TO NEW_NODE FOR THAT TREE.
!
! The fact that pointers are created undefined forces the use of this routine,
      type(data), pointer  :: tree
      type(state), pointer :: state_of_tree
!
      nullify(tree)
!
! but initialize some global variables too,
      count = 0
      last_code = 0
      current => tree
      parent => tree
!
! and set up state variable of tree.
      allocate(state_of_tree)
      state_of_tree%count     = count
      state_of_tree%last_code = last_code
      state_of_tree%current   => tree
      state_of_tree%parent    => parent
      tree_state              => state_of_tree
   end subroutine start

   subroutine update_state(tree)
! To update state variable in case next call to module is for a different
! tree.
      type(data), target, intent(in)  :: tree
      type(state), pointer:: state_of_tree
!
      state_of_tree           => tree%own_state
      state_of_tree%count     = count
      state_of_tree%last_code = last_code
      state_of_tree%current   => current
      state_of_tree%parent    => parent
   end subroutine update_state

   subroutine num_leaves(tree,nl)
      type(data), intent(in) :: tree
      integer nl

      nl=0
      call nl_calc(tree,nl)

   contains
      recursive subroutine nl_calc(tree,nl)
         type(data), intent(in) :: tree
         integer loop, nl

         ! Loop through whole tree
         if(size(tree%p) > 0) then
           do loop = 1, size(tree%p)
             if(associated(tree%p(loop)%pp)) call nl_calc (tree%p(loop)%pp,nl)
           end do
         else
           nl=nl+1
         endif

      end subroutine nl_calc

   end subroutine num_leaves

   subroutine rng_tree (tree,rngs)
   ! Traverse a complete tree or subtree, generating random numbers for edges probas
      use Dirichlet
      type(data), intent(in) :: tree
      integer ileaf
      real proba_parent
      real, dimension(:)     :: rngs

      ileaf=0
      proba_parent=1.0
      call rng_calc(tree,proba_parent,ileaf,rngs)

   contains
      recursive subroutine rng_calc(tree,proba_parent,ileaf,rngs)
         type(data), intent(in) :: tree
         integer            :: loop, nleafs, ileaf
         real               :: proba_parent
         real, allocatable  :: alpha(:), beta(:)
         real, dimension(:) :: rngs

         nleafs=size(tree%link)
         select case (tree%header)
           case('Drat')
             tree%user(1)%f_field=Drat(tree%y(1),tree%y(2))

           case('Disc')
             allocate(alpha(nleafs))
             alpha=tree%y(1:nleafs)
             tree%user(1)%f_field=Disc(alpha)
             deallocate(alpha)

           case('Diri')
             allocate(alpha(nleafs))
             alpha=tree%y(1:nleafs)
             allocate(beta(nleafs)); beta=0.0
             beta=tree%y(nleafs+1:2*nleafs)
             tree%user(1)%f_field=Diri(alpha,beta)
             deallocate(alpha,beta)

           case('Dirg')
             allocate(alpha(nleafs))
             alpha=tree%y(1:nleafs)
             allocate(beta(nleafs))
             beta=tree%y(nleafs+1:2*nleafs)
             tree%user(1)%f_field=Dirg(alpha,beta)
             deallocate(alpha,beta)

           case('Dirw')
             allocate(alpha(nleafs))
             alpha=tree%y(1:nleafs)
             allocate(beta(nleafs))
             beta=tree%y(nleafs+1:2*nleafs)
             tree%user(1)%f_field=Dirw(alpha,beta)
             deallocate(alpha,beta)

           case('Dirz')
             allocate(alpha(nleafs))
             alpha=tree%y(1:nleafs)
             allocate(beta(nleafs))
             beta=tree%y(nleafs+1:2*nleafs)
             tree%user(1)%f_field=Dirz(alpha,beta)
             deallocate(alpha,beta)

           case('Diut')
             allocate(alpha(nleafs))
             alpha=tree%y(1:nleafs)
             allocate(beta(nleafs))
             beta=tree%y(nleafs+1:2*nleafs)
             tree%user(1)%f_field=Diut(alpha,beta)
             deallocate(alpha,beta)

           case('Diun')
             tree%user(1)%f_field=Diun(nleafs)

           case('Dior')
             allocate(alpha(nleafs))
             alpha=tree%y(1:nleafs)
             tree%user(1)%f_field=Dior(alpha)
             deallocate(alpha)

          case('Logn')
             allocate(alpha(2))
             alpha=tree%y(1:2)
             tree%user(1)%f_field=Logn(alpha)
             deallocate(alpha)

           case('Logu')
             allocate(alpha(2))
             alpha=tree%y(1:2)
             tree%user(1)%f_field=Logu(alpha)
             deallocate(alpha)

          case('Norm')
             allocate(alpha(2))
             alpha=tree%y(1:2)
             tree%user(1)%f_field=Norm(alpha)
             deallocate(alpha)

           case('Nort')
             allocate(alpha(2))
             alpha=tree%y(1:2)
             allocate(beta(2))
             beta=tree%y(3:4)
             tree%user(1)%f_field=Nort(alpha,beta)
             deallocate(alpha,beta)

           case('Unif')
             allocate(alpha(2))
             alpha=tree%y(1:2)
             tree%user(1)%f_field=Unif(alpha)
             deallocate(alpha)

           case('Mlgn')
             allocate(alpha(nleafs))
             alpha=tree%y(1:nleafs)
             allocate(beta(nleafs))
             beta=tree%y(nleafs+1:2*nleafs)
             tree%user(1)%f_field= Mlogn(alpha,beta)
             deallocate(alpha,beta)

           case default
             tree%user(1)%f_field=1.

         end select
         ! Combine with parent's proba
         tree%user(1)%f_field = tree%user(1)%f_field * proba_parent

         ! Loop through whole tree
         if(size(tree%p)>0) then
           do loop = 1, size(tree%p)
             if(associated(tree%p(loop)%pp)) then
               proba_parent= tree%user(1)%f_field(loop)
               call rng_calc (tree%p(loop)%pp,proba_parent,ileaf,rngs)
             end if
           end do
         else
           ! Terminal leaf: store proba
           ileaf=ileaf+1
           rngs(ileaf)= proba_parent
         end if

      end subroutine rng_calc

   end subroutine rng_tree

   subroutine strip_spaces (chain)
   ! Strip spaces in chain
     IMPLICIT NONE
     CHARACTER (LEN=*), INTENT(inout) :: chain
     CHARACTER (LEN=len(chain))       :: str
     integer                          :: lc, i

     str=chain
     lc=len_trim(str)
     if(lc.le.1) return
     i=1
     !Heading spaces
     do while(str(i:i).eq.' ')
       str=str(i+1:lc)
       lc=lc-1
     enddo

     !Body spaces
     do while(i.lt.lc)
        if(str(i:i).eq.' ') then
           str=str(:i-1)//str(i+1:lc)
           lc=lc-1
        else
           i=i+1
        endif
     enddo
     chain=str

   end subroutine strip_spaces

   subroutine build_tree (chain,tree,nl)
   ! Parse the distribution expression and build a tree
     character (len=*)               :: chain
     type(data), intent(inout), pointer :: tree
     type(user_type)                 :: ran_tab(1)
     integer, intent(out)            :: nl
     integer                         :: i, nb_nodes, lc
     logical full_atomized
     type my_node
       character(len=chain_size)  :: chain
       character(len=4)           :: oper
       integer                    :: nlinks,links(max_leaves)
       real                       :: params(max_leaves)
     end type my_node
     type(my_node)                :: nodes(max_nodes)

     ! Init tree
     call start(tree)
     allocate(ran_tab(1)%f_field(max_leaves), &
              ran_tab(1)%i_field(max_leaves) )
     ran_tab(1)%f_field=0.

     call strip_spaces(chain)
     lc = len_trim(chain)
     if(lc.gt.chain_size) stop 'Pb dimensions: chain too small'

     nb_nodes=1
     nodes(nb_nodes)%chain = chain
     full_atomized=.FALSE.
     do while (.NOT.full_atomized)
        full_atomized=.TRUE.
        do i= 1, nb_nodes
           if(nodes(i)%chain .NE. '') then
              full_atomized=.FALSE.
              call atomize(i,nodes,nb_nodes)
              exit
           endif
        enddo
     enddo

     do i=  1, nb_nodes
       call new_node( tree, &
                      i, &
                      nodes(i)%oper, &
                      real_data=nodes(i)%params, &
                      links= (/ (nodes(i)%links(1:nodes(i)%nlinks)) /),&
                      user_data=ran_tab )
     enddo

     call num_leaves(tree,nl)
     print *,nl

     return

   contains

     subroutine atomize (indx,nodes,nb_nodes)
     ! Splits the distribution expression into elementary contributions

       character (len=4)               :: oper
       character (len=chain_size)      :: chain
       integer                         :: lc, lpo, lpc, indx, nel, i, j
       integer                         :: le,  nb_nodes, narg
       CHARACTER (len=chain_size)      :: elts(max_leaves)
       real, allocatable               :: alpha(:)
       type(my_node)                   :: nodes(:)

       !print *,'Atomizing node ',indx
       nel=0
       chain=nodes(indx)%chain
       lpo=index(chain,'(')
       lpc=len_trim(chain)
       oper=chain(1:4)
       chain=chain(lpo+1:lpc-1)
       lc = len_trim(chain)

       nodes(indx)%oper  = oper
       nodes(indx)%chain = ''

       select case (oper)

          case ('Dirg','Diut','Dirw','Dirz','Mlgn') ! Look for 2*nleafs params
            call get_arguments(chain,elts,nel)
            allocate(alpha(nel))
            do i= 1, nel
               le=len_trim(elts(i))
               if(index(elts(i),'*').gt.0) le=index(elts(i),'*')-1
               read(elts(i)(1:le),*) alpha(i)
               elts(i)=elts(i)(le+2:)
            enddo
            nel=nel/2 ! Upper elements are pure params
            nodes(indx)%nlinks= nel
            nodes(indx)%links = (/ (nb_nodes+j,j=1,nel) /)
            nodes(indx)%params= alpha
            deallocate(alpha)

          case ('Diri')        ! Look for nleafs+1/2 params
            call get_arguments(chain,elts,nel)
            allocate(alpha(nel))
            do i= 1, nel
               le=len_trim(elts(i))
               if(index(elts(i),'*').gt.0) le=index(elts(i),'*')-1
               read(elts(i)(1:le),*) alpha(i)
               elts(i)=elts(i)(le+2:)
            enddo
            nel=nel-1 ! Upper elements are pure params
            nodes(indx)%nlinks= nel
            nodes(indx)%links = (/ (nb_nodes+j,j=1,nel) /)
            nodes(indx)%params= alpha
            deallocate(alpha)

          case ('Diun')   ! Look for a single param
            read(chain,*) nel
            elts=''
            nodes(indx)%nlinks= nel
            nodes(indx)%links = (/ (nb_nodes+j,j=1,nel) /)

          case ('Dior')   ! Look for a single param or a vector
            call get_arguments(chain,elts,narg)

            if (narg==1) then
              ! Single argument : nb of leaves. Order is fixed.
              read(elts(1),*) nel
              nodes(indx)%nlinks= nel
              nodes(indx)%links = (/ (nb_nodes+j,j=1,nel) /)
              nodes(indx)%params= (/ (real(j)   ,j=1,nel) /)
              do i= 1, nel
                 elts(i)=''
              enddo

            else
              nel=narg
              allocate(alpha(nel))
              do i= 1, nel
                 le=len_trim(elts(i))
                 if(index(elts(i),'*').gt.0) le=index(elts(i),'*')-1
                 read(elts(i)(1:le),*) alpha(i)
                 elts(i)=elts(i)(le+2:)
              enddo
              nodes(indx)%nlinks= nel
              nodes(indx)%links = (/ (nb_nodes+j,j=1,nel) /)
              nodes(indx)%params= alpha
              deallocate(alpha)

            endif

          case ('Drat') ! Look for 2 params
            call get_arguments(chain,elts,nel)
            allocate(alpha(nel))
            do i= 1, nel
               le=len_trim(elts(i))
               if(index(elts(i),'*').gt.0) le=index(elts(i),'*')-1
               read(elts(i)(1:le),*) alpha(i)
               elts(i)=elts(i)(le+2:)
            enddo
            nodes(indx)%nlinks= nel
            nodes(indx)%links = (/ (nb_nodes+j,j=1,nel) /)
            nodes(indx)%params= alpha
            deallocate(alpha)

          case ('Disc')          ! Look for nleafs params
            call get_arguments(chain,elts,nel)
            allocate(alpha(nel))
            do i= 1, nel
               le=len_trim(elts(i))
               if(index(elts(i),'*').gt.0) le=index(elts(i),'*')-1
               read(elts(i)(1:le),*) alpha(i)
               elts(i)=elts(i)(le+2:)
            enddo
            nodes(indx)%nlinks= nel
            nodes(indx)%links = (/ (nb_nodes+j,j=1,nel) /)
            nodes(indx)%params= alpha
            deallocate(alpha)

          case ('Logn','Logu','Norm','Nort','Unif')
            call get_arguments(chain,elts,nel)
            allocate(alpha(nel))
            do i= 1, nel
               le=len_trim(elts(i))
               read(elts(i)(1:le),*) alpha(i)
            enddo
            nel=1
            elts=''
            nodes(indx)%nlinks= nel
            nodes(indx)%links = (/ (nb_nodes+j,j=1,nel) /)
            nodes(indx)%params= alpha
            deallocate(alpha)

       end select

       if(nb_nodes+nel.gt.max_nodes) stop 'Pb dimensions: max_nodes too small'
       do i= 1, nel
         oper='End'
         if(len_trim(elts(i)).gt.0) read(elts(i),'(A4)') oper
         nodes(nb_nodes+i)%oper  = oper
         if(oper .EQ. 'End') then
           nodes(nb_nodes+i)%chain = ''
         else
           nodes(nb_nodes+i)%chain = elts(i)
         endif
       enddo
       nb_nodes=nb_nodes + nel

     end subroutine atomize

     subroutine get_arguments(chain,elts,nel)
     ! Get the arguments of a distribution
     character (len=*)          :: chain
     integer                    :: j, npo, nel, j0, lv, lp, lc, ls
     CHARACTER (LEN=len(chain)) :: elts(max_leaves)

     nel=0
     j0=1
     lc = len_trim(chain)
     do while (j0.le.lc)
        lv=index(chain(j0:lc),',')
        ls=index(chain(j0:lc),';')
        if(ls*lv.gt.0) then
           lv=min(ls,lv)
        else
           lv=max(ls,lv)
        endif
        if(lv.ne.0) lv=lv+j0-1
        lp=index(chain(j0:lc),'(')
        if(lp.ne.0) lp=lp+j0-1
        if(lv.le.0) then ! Last element
          nel=nel+1
          if(nel.gt.max_leaves) stop 'Pb dimensions: max_leaves too small'
          elts(nel)=chain(j0:lc)
          exit
        else if(lp.le.0 .OR. lp.gt.lv) then ! Simple element
          nel=nel+1
          if(nel.gt.max_leaves) stop 'Pb dimensions: max_leaves too small'
          elts(nel)=chain(j0:lv-1)
          j0=lv+1
        else
          npo=1
          j=lp
          do while(npo.ne.0)
             j=j+1
             if(chain(j:j) == '(') then
               npo=npo+1
             else if(chain(j:j) == ')') then
               npo=npo-1
             endif
          enddo
          nel=nel+1
          if(nel.gt.max_leaves) stop 'Pb dimensions: max_leaves too small'
          elts(nel)=chain(j0:j)
          j0=j+2
        endif
     enddo

     end subroutine get_arguments

   end subroutine build_tree

   subroutine nestdiri(ndraw,chain,nl,sample)
   ! Generate sample from distribution expression
     use ddl
     implicit none
     integer                   :: i, nl, ndraw
     character(len=chain_size) :: chain
     real, allocatable         :: rngs(:)
     real                      :: sample(ndraw,max_leaves)
     type(data), pointer       :: a

     ! Initialize and build tree from chain
     nullify(a)
     call build_tree(chain,a,nl)
     allocate(rngs(nl))

     ! Generate sample
     do i= 1, ndraw
       call rng_tree(a,rngs)
       sample(i,1:nl)=rngs
     enddo

     ! Free space
     call finish(a)
     deallocate(rngs)

   end subroutine nestdiri

end module eagle

program test
  use ddl
  use eagle
  use Dirichlet
  implicit none
  integer                   :: i, iseed, is2, is1, nl, ndraw
  character(len=chain_size) :: chain
  real, allocatable         :: sample(:,:)

  ! Initialize RNG
  call system_clock(iseed)
  iseed=-iseed
  write(chain,*) iseed
  call phrtsd(chain,is1,is2)
  call setall(is1,is2)

  read(5,'(A)') chain
  read(chain,*) ndraw
  chain=chain(index(chain,' ')+1:)
  chain=trim(adjustl(chain)) ! remove spaces

  ! Generate sample
  allocate(sample(ndraw,max_leaves))
  call nestdiri(ndraw,chain,nl,sample)

  do i= 1, ndraw
    write(*,'(150e15.7)') sample(i,1:nl)
  enddo

  deallocate(sample)

end
