!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module global_variables
!-------10--------20--------30--------40--------50--------60--------70--------80

  integer,parameter :: Max_orbit=15
! constant
  real(8),parameter :: pi=3.1415926535897932d0
  real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0
  real(8) gaunt(16,16,16),nYlm(0:3,-3:3)
    

! log-mash parameter
  integer Nx
  real(8) Dx,Rmax,Rp
  real(8),allocatable :: xL(:),rL(:),expXL(:)


! DFT parameter 
  character(5) :: Functional
  real(8) :: rho_cut
! Perdew-Zunger LDA parameter 
  real(8),parameter :: gammaU=-0.1423d0,beta1U=1.0529d0
  real(8),parameter :: beta2U=0.3334d0,AU=0.0311d0,BU=-0.048d0
!  real(8),parameter :: CU=0.002d0,DU=-0.0116d0
  real(8),parameter :: CU=0.2019151940622859d-2, DU=-0.1163206637891297d-1 ! modified

! GS
  real(8),allocatable :: upsi(:,:),Vh(:),Eexc(:),Vexc(:),Veff(:)
  real(8),allocatable :: VL(:),Vnucl(:)
  real(8),allocatable :: rho(:)
  real(8),allocatable :: occA(:),esp(:)
  integer Nupsi
  integer,allocatable :: L_upsi(:)
  real(8) :: total_energy
  real(8),parameter :: rate=0.1d0
  integer :: Nscf
! atom date
  integer ZA
  character(5) atom_name
  real(8),allocatable :: occ_input(:)
  integer,allocatable :: orbit_L(:)
  character(10),allocatable :: orbit_state(:)
  character(2),allocatable :: orbit_name(:),name_upsi(:)
  
! PAW
  integer NuPR
  integer Nunocc(0:3)
  integer,allocatable :: L_proj(:)
  real(8),allocatable :: Cr_phi_P(:,:)
  real(8),allocatable :: Rcut(:)
  integer,allocatable :: iRc(:)
  real(8) Rfilt,Rcomp,Rcore,RcutPR
  integer iRcf,iRcomp,iRcore,iRcPR
  real(8),allocatable :: uAE(:,:),uPS(:,:),uPR(:,:)
  real(8),allocatable :: nc_A(:),nc_P(:),V_bar(:),gL(:,:)
    

end module global_variables
!-------10--------20--------30--------40--------50--------60--------70--------80
program main
  use global_variables
  implicit none
  integer iter,p,ix
  integer i,j

! calc  Ground state ------------------------------------------------!
  call prep_calc_parameter
  call init_wf
  call wf_normalize
  call upsi_rho
  call Hartree(0)
  call Exc_Cor_PZ_LDA

  do iter=1,400
    call solve_Kohn_Sham
    call wf_normalize
    call upsi_rho
    call Hartree(0)
    call Exc_Cor_PZ_LDA
  end do

  call calc_total_energy
  write(*,"(A,2x,f12.4)")'Total energy (a.u.)',total_energy
  do p=1,Nupsi
    write(*,"(I5,2x,A,2x,f14.6,2x,f16.4)")p,name_upsi(p),esp(p),esp(p)*2d0*Ry
  end do

  write(*,*)'End prepare'
  write(*,*)'Rmax=',Rmax
  write(*,*)'Dx=',Dx,1d0/32d0
  write(*,*)'Nx=',Nx
  write(*,*)'Rp=',Rp
  
  call write_wavefunction
    
! end calc  Ground state --------------------------------------------!
  
end program main
!-------10--------20--------30--------40--------50--------60----------72
subroutine prep_calc_parameter
  use global_variables
  implicit none
  integer i,j,ix,m
  
  call prep_orbit_L
  
  allocate(occ_input(Max_orbit),orbit_state(Max_orbit))
    
  
! input file

  read(*,*)atom_name   ! atom name
  read(*,*)ZA          ! nuclear chage
  read(*,*)Nunocc(:)   !number unoccupied orbit
  do i=1,Max_orbit
    read(*,*)occ_input(i),orbit_state(i) ! occupation & state
  end do
  read(*,*)Dx           ! grid size
  read(*,*)Nx          ! mesh points
  read(*,*)Rmax        ! radious simulation sphere
  read(*,*)rho_cut        ! rho_cut
! count the number
  NuPR=0
  Nupsi=0
  do i=1,Max_orbit
    if((orbit_state(i)=='v')) then
      NuPR=NuPR+1
      Nupsi=Nupsi+1
    end if
    if(orbit_state(i)=='c')Nupsi=Nupsi+1
  end do
  NuPR=NuPR+sum(Nunocc)
  Nupsi=Nupsi+sum(Nunocc)
  
  Rp=Rmax/(exp(Dx*dble(Nx))-1d0)
  
  allocate(xL(0:Nx),rL(0:Nx),expXL(0:Nx))
  do ix=0,Nx
    xL(ix)=Dx*dble(ix)
    rL(ix)=Rp*(exp(Dx*dble(ix))-1d0)
    expXL(ix)=exp(Dx*dble(ix))
  end do
  
    
  allocate(L_proj(NuPR),L_upsi(Nupsi),name_upsi(Nupsi))
  allocate(occA(Nupsi),esp(Nupsi))
  
  NuPR=0
  Nupsi=0
  do i=1,Max_orbit
    if((orbit_state(i)=='v'))then
      NuPR=NuPR+1
      Nupsi=Nupsi+1
      L_proj(NuPR)=orbit_L(i)
      L_upsi(Nupsi)=orbit_L(i)
      name_upsi(Nupsi)=orbit_name(i)
      occA(Nupsi)=occ_input(i)
    end if
    if(orbit_state(i)=='c')then
      Nupsi=Nupsi+1
      L_upsi(Nupsi)=orbit_L(i)
      name_upsi(Nupsi)=orbit_name(i)
      occA(Nupsi)=occ_input(i)
    end if
  end do
  
  do i=0,3
    do j=1,Nunocc(i)
      NuPR=NuPR+1
      Nupsi=Nupsi+1
      L_proj(NuPR)=i
      L_upsi(Nupsi)=i
      select case(i)
      case(0)
         name_upsi(Nupsi)="us"
      case(1)
         name_upsi(Nupsi)="up"
      case(2)
         name_upsi(Nupsi)="ud"
      case(3)
         name_upsi(Nupsi)="uf"
      end select
      occA(Nupsi)=0d0
    end do
  end do

  allocate(upsi(0:Nx,Nupsi),rho(0:Nx))
  rho=0d0
  allocate(Vh(0:Nx),Vexc(0:Nx),Veff(0:Nx),Eexc(0:Nx))
  allocate(VL(0:Nx),Vnucl(0:Nx))
  
  Vnucl(1:Nx)=-dble(ZA)/rL(1:Nx);Vnucl(0)=0d0
  VL(1:Nx)=0.5d0/(rL(1:Nx)**2);VL(0)=0d0
  
  return
end subroutine prep_calc_parameter
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine prep_orbit_L
  use global_variables
  implicit none
  
  allocate(orbit_L(Max_orbit),orbit_name(Max_orbit))

  orbit_L(1)= 0 ;orbit_name(1)= "1s"
  orbit_L(2)= 0 ;orbit_name(2)= "2s"
  orbit_L(3)= 1 ;orbit_name(3)= "2p"
  orbit_L(4)= 0 ;orbit_name(4)= "3s"
  orbit_L(5)= 1 ;orbit_name(5)= "3p"
  orbit_L(6)= 2 ;orbit_name(6)= "3d"
  orbit_L(7)= 0 ;orbit_name(7)= "4s"
  orbit_L(8)= 1 ;orbit_name(8)= "4p"
  orbit_L(9)= 2 ;orbit_name(9)= "4d"
  orbit_L(10)=3 ;orbit_name(10)="4f"
  orbit_L(11)=0 ;orbit_name(11)="5s"
  orbit_L(12)=1 ;orbit_name(12)="6s"
  orbit_L(13)=2 ;orbit_name(13)="5d"
  orbit_L(14)=0 ;orbit_name(14)="6s"
  orbit_L(15)=1 ;orbit_name(15)="6p"
  return
end subroutine prep_orbit_L
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine init_wf
  use global_variables
  implicit none
  real(8) r
  integer ix,p
  
  do p=1,Nupsi
    upsi(:,p)=rL(:)**(L_upsi(p)+1)*exp(-5.125*rL(:)**2)
  end do
  
  return
end subroutine init_wf
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine wf_normalize
  use global_variables
  implicit none
  integer p
  real(8) S
  
  do p=1,Nupsi
    S=sum(upsi(:,p)**2*expXL(:))*Dx*Rp
    S=1d0/sqrt(S)
    upsi(:,p)=S*upsi(:,p)
  end do
  
  return
end subroutine wf_normalize
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine upsi_rho
  use global_variables
  implicit none
  real(8),allocatable :: rho_old(:)
  integer p
  
  allocate(rho_old(0:Nx))
  rho_old=rho
  rho=0d0
  do p=1,Nupsi
    rho(1:Nx)=rho(1:Nx)+occA(p)/(4d0*pi)*(upsi(1:Nx,p)/rL(1:Nx))**2
  end do
  rho(0)=2d0*rho(1)-rho(2)
  rho=rho*rate+(1d0-rate)*rho_old
  return
end subroutine upsi_rho
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine Hartree(l)
  use global_variables
  implicit none
  integer ix,l
  real(8),allocatable :: Sn(:),Vn(:),Kn(:)
  real(8) Qrho,Const
  real(8) fact1,fact2,fact3
  real(8) temp1,temp2
  
  allocate(Sn(0:Nx),Vn(0:Nx),Kn(0:Nx))
  
  Qrho=sum(expXL*rL**(2+l)*rho)*Rp*Dx*4d0*pi/dble(2*l+1)
  Sn=-4d0*pi*Rp**2*rL*exp(1.5d0*xL)*rho
  fact1=Dx**2/12d0
  fact2=5d0*Dx**2/12d0
  Vn(0:2)=rL(0:2)**(l+1)*exp(-0.5d0*xL(0:2))
  do ix=1,Nx
    Kn(ix)=-0.25d0-(Rp*expXL(ix)/rL(ix))**2*dble(l*(l+1))
  end do
  
  do ix=3,Nx
    temp1=2d0*(1d0-fact2*Kn(ix-1))*Vn(ix-1) &
      -(1d0+fact1*Kn(ix-2))*Vn(ix-2)
    temp2=fact1*(Sn(ix)+10d0*Sn(ix-1)+Sn(ix-2))
    Vn(ix)=(temp1+temp2)/(1d0+fact1*Kn(ix))
  end do
  
  Vh(1:Nx)=Vn(1:Nx)*exp(0.5d0*xL(1:Nx))/rL(1:Nx)
  Const=Qrho/(rL(Nx)**(2*l+1))-Vh(Nx)/(rL(Nx)**l)
  Vh(1:Nx)=Vh(1:Nx)+Const*rL(1:Nx)**l
  Vh(0)=2d0*Vh(1)-Vh(2)
  return
end subroutine Hartree
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine Exc_Cor_PZ_LDA
  use global_variables
  implicit none
  integer ix
  real(8) trho,rs,rssq,rsln,V_xc,E_xc
  real(8),parameter :: const = 0.75d0*(3d0/(2d0*pi))**(2d0/3d0)
  
  do ix=0,Nx
    trho=rho(ix)+1d-20
    rs=(3d0/(4*Pi*trho))**(1d0/3d0)
    V_xc=-4d0/3d0*const/rs
    E_xc=-const/rs
    if(rs>1d0) then
      rssq=sqrt(rs)
      V_xc=V_xc+gammaU*(1d0+7d0/6d0*beta1U*rssq+4d0/3d0*beta2U*rs) &
        /(1d0+beta1U*rssq+beta2U*rs)**2
      E_xc=E_xc+gammaU/(1d0+beta1U*rssq+beta2U*rs)
    else
      rsln=log(rs)
      V_xc=V_xc+AU*rsln+(BU-AU/3d0) &
        +2d0/3d0*CU*rs*rsln+(2d0*DU-CU)/3d0*rs
      E_xc=E_xc+AU*rsln+BU+CU*rs*rsln+DU*rs
    endif
    Vexc(ix)=V_xc
    Eexc(ix)=E_xc
  end do
  return
end subroutine Exc_Cor_PZ_LDA
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine solve_Kohn_Sham
  use global_variables
  implicit none
  integer Node_L(0:3),node,p,ix,n
  real(8) epsilon,ref_epsilon,epsilon_min,epsilon_max
  real(8) fact1,temp1,temp2
  real(8),allocatable :: v_upsi(:,:),K2(:)
  allocate(v_upsi(0:Nx,Nupsi),K2(0:Nx))

  Node_L=0
  ref_epsilon=1d-12
  fact1=Dx**2/12d0
  
  do p=1,Nupsi
    v_upsi(0:2,p)=rL(0:2)**(L_upsi(p)+1)*exp(-0.5d0*XL(0:2))
    Veff=Vh+Vexc+Vnucl+VL*dble(L_upsi(p)*(L_upsi(p)+1))
    epsilon_min=-0.5d0*dble(ZA)**2-10d0
    epsilon_max=3d0       
! search eigen value -------------------------------------------------!
    do
      epsilon=0.5d0*(epsilon_max+epsilon_min)
      K2=-0.25d0+2d0*((Rp*expXL)**2)*(epsilon-Veff)
      node=0
      do ix=3,Nx
        temp1=2d0*(1d0-5d0*fact1*K2(ix-1))*v_upsi(ix-1,p) &
          -(1d0+fact1*K2(ix-2))*v_upsi(ix-2,p)
        temp2=(1d0+fact1*K2(ix))
        v_upsi(ix,p)=temp1/temp2
        if(v_upsi(ix-1,p)*v_upsi(ix,p)<0d0)node=node+1
        if(abs(v_upsi(ix,p))>=1d2)exit
      end do
      
      if(node>node_L(L_upsi(p)))then
        epsilon_max=epsilon
      else
        epsilon_min=epsilon
      end if
      
      if(epsilon_max-epsilon_min<=ref_epsilon)exit
    end do
! end search eigen value ---------------------------------------------!
! calc wave function -------------------------------------------------!
    epsilon=epsilon_max
    esp(p)=epsilon
    K2=-0.25+2d0*(Rp*expXL)**2*(epsilon-Veff)
    node=0
    do ix=3,Nx
      temp1=2d0*(1d0-5d0*fact1*K2(ix-1))*v_upsi(ix-1,p) &
        -(1d0+fact1*K2(ix-2))*v_upsi(ix-2,p)
      temp2=(1d0+fact1*K2(ix))
      v_upsi(ix,p)=temp1/temp2
      if(v_upsi(ix,p)*v_upsi(ix-1,p)<0d0)node=node+1
      if(node>node_L(L_upsi(p)))exit
    end do
    v_upsi(ix:Nx,p)=0d0
! end calc wave function ---------------------------------------------!
    node_L(L_upsi(p))=node_L(L_upsi(p))+1
  end do
  
  do p=1,Nupsi
    upsi(:,p)=v_upsi(:,p)*exp(0.5d0*XL(:))
  end do
  return
end subroutine solve_Kohn_Sham
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_total_energy
  use global_variables
  implicit none
  integer p

  total_energy = 0d0
  do p = 1,Nupsi
     total_energy = total_energy &
          + occA(p)*esp(p)
  end do

  total_energy = total_energy + 4d0*pi*sum(rL(:)**2*rho(:)*expXL(:)*(Eexc(:)-Vexc(:)-0.5*Vh(:)))*Dx*Rp

  return
end subroutine calc_total_energy
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine write_wavefunction
  use global_variables
  implicit none
  integer :: ix

  open(103,file='results.dat')
  do ix=0,Nx
    write(103,'(999e26.16e3)')rL(ix),rho(ix),Vnucl(ix),Vh(ix),Vexc(ix) &
         ,upsi(ix,:)
  end do
  close(103)

  return
end subroutine write_wavefunction

!-------10--------20--------30--------40--------50--------60--------70--------80
