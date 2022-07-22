      program testgam
      character*80 phrase
      call init_io
      call get_random_seed
      read (*,*) x
      do i= 1, 1000
         print *,i,gengam(1.0,x)
      enddo
      end
