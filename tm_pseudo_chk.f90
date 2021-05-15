module global_variables
  real(8),parameter :: pi=3.1415926535897932d0


  ! pseudotepontial parameters
  real(8) :: Znucl, r_c
  ! coefficient for pseudopotentials
  real(8) :: cpp0,cpp2,cpp4,cpp6,cpp8,cpp10,cpp12
  real(8) :: cpp2_max, cpp2_min, norm_diff_max, norm_diff_min
  real(8) :: cpp0_max, cpp0_min

  real(8) :: Amat_pp(5,5),bvec(5),cvec(5),Bmat_pp(5,5)

  interface matrix_inverse
     module procedure matrix_inverse_double
!     module procedure matrix_inverse_complex
  end interface matrix_inverse

  contains

subroutine matrix_inverse_double(rmat)
  implicit none
  real(8),intent(inout) :: rmat(:,:)
  integer :: nn
! for lapack
  integer :: lwork, info
  integer, allocatable :: ipiv(:) ! dimension N
  real(8), allocatable :: work(:) ! dimension LWORK  


  nn = size(rmat,dim=1)
  lwork = nn * max(nn, 64)
    
  allocate(ipiv(nn),work(lwork))

  call dgetrf(nn, nn, rmat, nn, ipiv, info)  ! factorize
  call dgetri(nn, rmat, nn, ipiv, work, lwork, info)  ! inverse
  
  deallocate(ipiv,work)

end subroutine matrix_inverse_double
  
end module global_variables
!-----------------------------------------------------------------!
program main
  use global_variables
  implicit none

  ! ps parameteres
  Znucl = 1d0
  r_c = 0.1d0

  !  call parameter_search
    call parameter_search_TM_mod
  call write_output
  
  
end program main
!-----------------------------------------------------------------!
subroutine parameter_search
  use global_variables
  implicit none
  real(8) :: norm

  Amat_pp = 0d0
  Amat_pp(1,1) = 1d0; Amat_pp(1,2) = r_c**6; Amat_pp(1,3) = r_c**8
  Amat_pp(1,4) = r_c**10; Amat_pp(1,5) = r_c**12

  Amat_pp(2,1) = 0d0; Amat_pp(2,2) = 6d0*r_c**5; Amat_pp(2,3) = 8d0*r_c**7
  Amat_pp(2,4) = 10d0*r_c**9; Amat_pp(2,5) = 12d0*r_c**11

  Amat_pp(3,1) = 0d0; Amat_pp(3,2) = 30d0*r_c**4; Amat_pp(3,3) = 56d0*r_c**6
  Amat_pp(3,4) = 90d0*r_c**8; Amat_pp(3,5) = 132d0*r_c**10
  
  Amat_pp(4,1) = 0d0; Amat_pp(4,2) = 120d0*r_c**3; Amat_pp(4,3) = 336d0*r_c**5
  Amat_pp(4,4) = 720d0*r_c**7; Amat_pp(4,5) = 1320d0*r_c**9

  Amat_pp(5,1) = 0d0; Amat_pp(5,2) = 360d0*r_c**2; Amat_pp(5,3) = 1680d0*r_c**4
  Amat_pp(5,4) = 5040d0*r_c**6; Amat_pp(5,5) = 11880d0*r_c**8

  Bmat_pp = Amat_pp
  call matrix_inverse(Bmat_pp)

  norm_diff_max = -1d0
  norm_diff_min =  1d0
  cpp2 = 1d0
  do
     cpp4 = -cpp2/5d0
     bvec(1) = log(sqrt(4d0*Znucl**3)*exp(-znucl*r_c))-cpp2*r_c**2-cpp4*r_c**4
     bvec(2) = -znucl-2d0*cpp2*r_c-4d0*cpp4*r_c**3
     bvec(3) = -2d0*cpp2-12d0*cpp4*r_c**2
     bvec(4) = -24d0*cpp4*r_c
     bvec(5) = -24d0*cpp4
     cvec = matmul(Bmat_pp,bvec)
     cpp0 = cvec(1)
     cpp6 = cvec(2)
     cpp8 = cvec(3)
     cpp10 = cvec(4)
     cpp12 = cvec(5)

     call calc_norm(norm)
     if(norm > 0d0)then
        if(norm_diff_max <= 0d0)then
           norm_diff_max = norm
           cpp2_max = cpp2
        else if(norm_diff_max > norm)then
           norm_diff_max = norm
           cpp2_max = cpp2
        end if
     else if(norm < 0d0)then
        if(norm_diff_min >= 0d0)then
           norm_diff_min = norm
           cpp2_min = cpp2
        else if(norm_diff_min < norm)then
           norm_diff_min = norm
           cpp2_min = cpp2
        end if
     end if

     if(norm_diff_max > 0d0 .and. norm_diff_min < 0d0)exit
     cpp2 = -2d0*cpp2
     
  end do


  ! bi-section
  do
     cpp2 = 0.5d0*(cpp2_max + cpp2_min)
     cpp4 = -cpp2/5d0
     bvec(1) = log(sqrt(4d0*Znucl**3)*exp(-znucl*r_c))-cpp2*r_c**2-cpp4*r_c**4
     bvec(2) = -znucl-2d0*cpp2*r_c-4d0*cpp4*r_c**3
     bvec(3) = -2d0*cpp2-12d0*cpp4*r_c**2
     bvec(4) = -24d0*cpp4*r_c
     bvec(5) = -24d0*cpp4
     cvec = matmul(Bmat_pp,bvec)
     cpp0 = cvec(1)
     cpp6 = cvec(2)
     cpp8 = cvec(3)
     cpp10 = cvec(4)
     cpp12 = cvec(5)

     call calc_norm(norm)

     if(norm >= 0d0)then
        if(norm_diff_max >= norm)then
           norm_diff_max = norm
           cpp2_max = cpp2
        else
           write(*,*)"Error in bi-section +"
           stop
        end if
     else
        if(norm_diff_min <= norm)then
           norm_diff_min = norm
           cpp2_min = cpp2
        else
           write(*,*)"Error in bi-section -"
           stop
        end if        
     end if

     if(abs(cpp2_max-cpp2_min)<1d-12)exit
        
  end do


  cpp2 = 0.5d0*(cpp2_max + cpp2_min)
  cpp4 = -cpp2/5d0
  bvec(1) = log(sqrt(4d0*Znucl**3)*exp(-znucl*r_c))-cpp2*r_c**2-cpp4*r_c**4
  bvec(2) = -znucl-2d0*cpp2*r_c-4d0*cpp4*r_c**3
  bvec(3) = -2d0*cpp2-12d0*cpp4*r_c**2
  bvec(4) = -24d0*cpp4*r_c
  bvec(5) = -24d0*cpp4
  cvec = matmul(Bmat_pp,bvec)
  cpp0 = cvec(1)
  cpp6 = cvec(2)
  cpp8 = cvec(3)
  cpp10 = cvec(4)
  cpp12 = cvec(5)

  write(*,*)"cpp0=",cpp0
  write(*,*)"cpp2=",cpp2
  write(*,*)"cpp4=",cpp4
  write(*,*)"cpp6=",cpp6
  write(*,*)"cpp8=",cpp8
  write(*,*)"cpp10=",cpp10
  write(*,*)"cpp12=",cpp12

  call calc_norm(norm)
  write(*,*)"norm=",norm
  
end subroutine parameter_search
!-----------------------------------------------------------------!
subroutine parameter_search_TM_mod
  use global_variables
  implicit none
  real(8) :: norm
  real(8) :: Amat_pp_m(6,6),bvec_m(6),cvec_m(6),Bmat_pp_m(6,6)

  Amat_pp_m = 0d0
  Amat_pp_m(1,1) = r_c**2; Amat_pp_m(1,2) = r_c**4; Amat_pp_m(1,3) = r_c**6
  Amat_pp_m(1,4) = r_c**8; Amat_pp_m(1,5) = r_c**10; Amat_pp_m(1,6) = r_c**12

  Amat_pp_m(2,1) = 2*r_c; Amat_pp_m(2,2) = 4*r_c**3; Amat_pp_m(2,3) = 6*r_c**5
  Amat_pp_m(2,4) = 8*r_c**7; Amat_pp_m(2,5) = 10*r_c**9; Amat_pp_m(2,6) = 12*r_c**11

  Amat_pp_m(3,1) = 2; Amat_pp_m(3,2) = 4*3*r_c**2; Amat_pp_m(3,3) = 6*5*r_c**4
  Amat_pp_m(3,4) = 8*7*r_c**6; Amat_pp_m(3,5) = 10*9*r_c**8; Amat_pp_m(3,6) = 12*11*r_c**10

  Amat_pp_m(4,1) = 0d0; Amat_pp_m(4,2) = 4*3*2*r_c; Amat_pp_m(4,3) = 6*5*4*r_c**3
  Amat_pp_m(4,4) = 8*7*6*r_c**5; Amat_pp_m(4,5) = 10*9*8*r_c**7; Amat_pp_m(4,6) = 12*11*10*r_c**9

  Amat_pp_m(5,1) = 0d0; Amat_pp_m(5,2) = 4*3*2; Amat_pp_m(5,3) = 6*5*4*3*r_c**2
  Amat_pp_m(5,4) = 8*7*6*5*r_c**4; Amat_pp_m(5,5) = 10*9*8*7*r_c**6
  Amat_pp_m(5,6) = 12*11*10*9*r_c**8

  Amat_pp_m(6,1) = 0d0; Amat_pp_m(6,2) = 0d0; Amat_pp_m(6,3) = 6*5*4*3*2*r_c
  Amat_pp_m(6,4) = 8*7*6*5*4*r_c**3; Amat_pp_m(6,5) = 10*9*8*7*6*r_c**5
  Amat_pp_m(6,6) = 12*11*10*9*8*r_c**7  
  
  Bmat_pp_m = Amat_pp_m
  call matrix_inverse(Bmat_pp_m)

  norm_diff_max = -1d0
  norm_diff_min =  1d0
  cpp0 = 1d0
  do
     bvec_m(1) = log(sqrt(4d0*Znucl**3)*exp(-znucl*r_c))-cpp0
     bvec_m(2) = -znucl
     bvec_m(3) = 0d0
     bvec_m(4) = 0d0
     bvec_m(5) = 0d0
     bvec_m(6) = 0d0
     cvec_m = matmul(Bmat_pp_m,bvec_m)
     cpp2 = cvec_m(1)
     cpp4 = cvec_m(2)
     cpp6 = cvec_m(3)
     cpp8 = cvec_m(4)
     cpp10 = cvec_m(5)
     cpp12 = cvec_m(6)

     call calc_norm(norm)
     if(norm > 0d0)then
        if(norm_diff_max <= 0d0)then
           norm_diff_max = norm
           cpp0_max = cpp0
        else if(norm_diff_max > norm)then
           norm_diff_max = norm
           cpp0_max = cpp0
        end if
     else if(norm < 0d0)then
        if(norm_diff_min >= 0d0)then
           norm_diff_min = norm
           cpp0_min = cpp0
        else if(norm_diff_min < norm)then
           norm_diff_min = norm
           cpp0_min = cpp0
        end if
     end if

     if(norm_diff_max > 0d0 .and. norm_diff_min < 0d0)exit
     cpp0 = -2d0*cpp0
     
  end do


  ! bi-section
  do
     cpp0 = 0.5d0*(cpp0_max + cpp0_min)

     bvec_m(1) = log(sqrt(4d0*Znucl**3)*exp(-znucl*r_c))-cpp0
     bvec_m(2) = -znucl
     bvec_m(3) = 0d0
     bvec_m(4) = 0d0
     bvec_m(5) = 0d0
     bvec_m(6) = 0d0
     cvec_m = matmul(Bmat_pp_m,bvec_m)
     cpp2 = cvec_m(1)
     cpp4 = cvec_m(2)
     cpp6 = cvec_m(3)
     cpp8 = cvec_m(4)
     cpp10 = cvec_m(5)
     cpp12 = cvec_m(6)     

     call calc_norm(norm)

     if(norm >= 0d0)then
        if(norm_diff_max >= norm)then
           norm_diff_max = norm
           cpp0_max = cpp0
        else
           write(*,*)"Error in bi-section +"
           stop
        end if
     else
        if(norm_diff_min <= norm)then
           norm_diff_min = norm
           cpp0_min = cpp0
        else
           write(*,*)"Error in bi-section -"
           stop
        end if        
     end if

     if(abs(cpp0_max-cpp0_min)<1d-12)exit
        
  end do


  cpp0 = 0.5d0*(cpp0_max + cpp0_min)
  bvec_m(1) = log(sqrt(4d0*Znucl**3)*exp(-znucl*r_c))-cpp0
  bvec_m(2) = -znucl
  bvec_m(3) = 0d0
  bvec_m(4) = 0d0
  bvec_m(5) = 0d0
  bvec_m(6) = 0d0
  cvec_m = matmul(Bmat_pp_m,bvec_m)
  cpp2 = cvec_m(1)
  cpp4 = cvec_m(2)
  cpp6 = cvec_m(3)
  cpp8 = cvec_m(4)
  cpp10 = cvec_m(5)
  cpp12 = cvec_m(6)     
  
  call calc_norm(norm)  



  write(*,*)"cpp0=",cpp0
  write(*,*)"cpp2=",cpp2
  write(*,*)"cpp4=",cpp4
  write(*,*)"cpp6=",cpp6
  write(*,*)"cpp8=",cpp8
  write(*,*)"cpp10=",cpp10
  write(*,*)"cpp12=",cpp12

  call calc_norm(norm)
  write(*,*)"norm=",norm
  
end subroutine parameter_search_TM_mod
!-----------------------------------------------------------------!
subroutine calc_norm(norm)
  use global_variables
  implicit none
  real(8),intent(out)  :: norm
  integer,parameter :: Nr_norm = 10000
  real(8) :: rr, dr, R_PP, R_AE
  integer :: ir


  dr = r_c/Nr_norm
  
  norm = 0d0
  do ir = 1,Nr_norm-1
     rr = dr*ir
     R_AE = sqrt(4d0*Znucl**3)*exp(-Znucl*rr)
     R_PP = exp(cpp0 + cpp2*rr**2 + cpp4*rr**4 + cpp6*rr**6 &
              + cpp8*rr**8 + cpp10*rr**10 + cpp12*rr**12)
     
     norm = norm + dr*rr**2*(R_PP**2 - R_AE**2)
  end do
  rr = dr*Nr_norm
  R_AE = sqrt(4d0*Znucl**3)*exp(-Znucl*rr)
  R_PP = exp(cpp0 + cpp2*rr**2 + cpp4*rr**4 + cpp6*rr**6 &
           + cpp8*rr**8 + cpp10*rr**10 + cpp12*rr**12)
  
  norm = norm + 0.5d0*dr*rr**2*(R_PP**2 - R_AE**2)
  

end subroutine calc_norm
!-----------------------------------------------------------------!
subroutine write_output
  use global_variables
  implicit none
  integer,parameter :: Nr_out = 20000
  real(8),parameter :: rmax = 10d0, dr = rmax/Nr_out
  integer :: ir
  real(8) :: rr, R_PP, R_AE, V_PP, V_AE
  real(8) :: P1_r, P2, P1

  open(20,file="wf_pot.out")
  do ir = 0, Nr_out
     rr = dr*ir
     R_AE = sqrt(4d0*Znucl**3)*exp(-znucl*rr)
     if(rr<=r_c)then
        R_PP = exp(cpp0 + cpp2*rr**2 + cpp4*rr**4 + cpp6*rr**6 &
             + cpp8*rr**8 + cpp10*rr**10 + cpp12*rr**12)
        P1_r = 2d0*cpp2 + 4d0*cpp4*rr**2 + 6d0*cpp6*rr**4 + 8d0*cpp8*rr**6 &
              +10d0*cpp10*rr**8 + 12d0*cpp12*rr**10
        P1 = P1_r*rr
        P2 = 2d0*cpp2 + 12d0*cpp4*rr**2 + 30d0*cpp6*rr**4 + 56d0*cpp8*rr**6 &
             +90d0*cpp10*rr**8 + 132d0*cpp12*rr**10
        V_PP = -0.5d0*znucl**2 +P1_r +0.5d0*(P2 + P1**2)
        if(rr == 0d0)then
           V_AE = -znucl/(dr/10d0)
        else
           V_AE = -znucl/rr
        end if        
        
     else
        R_PP = R_AE
        V_AE = -znucl/rr
        V_PP = V_AE

     end if

     write(20,"(999e26.16e3)")rr,R_AE,R_PP,V_AE,V_PP
     
  end do
  close(20)
  

  
end subroutine write_output
!-----------------------------------------------------------------!
