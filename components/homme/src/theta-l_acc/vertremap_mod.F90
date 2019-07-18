#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module vertremap_mod
  use vertremap_base, only: remap1, remap1_nofilter,remap1_1_openacc,remap1_3_openacc

  use kinds, only                  : real_kind,int_kind
  use dimensions_mod, only         : np,nlev,qsize,nlevp,npsq,nelemd
  use hybvcoord_mod, only          : hvcoord_t
  use element_mod, only            : element_t
  use perf_mod, only               : t_startf, t_stopf  ! _EXTERNAL
  use parallel_mod, only           : abortmp, parallel_t
  use control_mod, only : vert_remap_q_alg
  use element_ops, only : set_theta_ref
  use eos, only : get_phinh,get_phinh_openacc
  implicit none
  private
  public :: vertical_remap

contains


  subroutine vertical_remap(hybrid,elem,hvcoord,dt,np1,np1_qdp,nets,nete)

  ! This routine is called at the end of the vertically Lagrangian
  ! dynamics step to compute the vertical flux needed to get back
  ! to reference eta levels
  !
  ! input:
  !     derived%dp()  delta p on levels at beginning of timestep
  !     state%dp3d(np1)  delta p on levels at end of timestep
  ! output:
  !     state%ps_v(np1)          surface pressure at time np1
  !     derived%eta_dot_dpdn()   vertical flux from final Lagrangian
  !                              levels to reference eta levels
  !
  use kinds,          only: real_kind
  use hybvcoord_mod,  only: hvcoord_t
  use control_mod,    only: rsplit
  use hybrid_mod,     only: hybrid_t
  use physical_constants, only : Cp,g

  type (hybrid_t),  intent(in)    :: hybrid  ! distributed parallel structure (shared)
  type (element_t), intent(inout) :: elem(:)
  type (hvcoord_t)                :: hvcoord
  real (kind=real_kind)           :: dt

  integer :: ie,i,j,k,np1,nets,nete,np1_qdp
  integer :: q
  
  real (kind=real_kind), dimension(nelemd,np,np,nlev)  :: dp_ie,dp_star_ie
  real (kind=real_kind), dimension(nelemd,np,np,nlevp) :: phi_ref_ie
  real (kind=real_kind), dimension(nelemd,np,np,nlev,5)  :: ttmp_ie
  real (kind=real_kind), dimension(np,np) :: tmp_arr

  call t_startf('vertical_remap')

  ! reference levels:
  !   dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v(i,j)
  !   hybi(1)=0          pure pressure at top of atmosphere
  !   hyai(1)=ptop
  !   hyai(nlev+1) = 0   pure sigma at bottom
  !   hybi(nlev+1) = 1
  !
  ! sum over k=1,nlev
  !  sum(dp(k)) = (hyai(nlev+1)-hyai(1))*ps0 + (hybi(nlev+1)-hybi(1))*ps_v
  !             = -ps0 + ps_v
  !  ps_v =  ps0+sum(dp(k))
  !
  ! reference levels:
  !    dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v
  ! floating levels:
  !    dp_star(k) = dp(k) + dt_q*(eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  ! hence:
  !    (dp_star(k)-dp(k))/dt_q = (eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  !
#ifdef OPENACC_HOMME
!$acc enter data create(dp_ie,dp_star_ie, phi_ref_ie, ttmp_ie)
!$acc parallel loop gang vector present(elem,hvcoord) private(tmp_arr)
#endif  
   do ie=nets,nete
     tmp_arr(:,:) = 0.0
     do k=1,nlev
        tmp_arr(:,:) = tmp_arr(:,:)+ elem(ie)%state%dp3d(:,:,k,np1)
     enddo
     ! update final ps_v
     elem(ie)%state%ps_v(:,:,np1) = hvcoord%hyai(1)*hvcoord%ps0 + tmp_arr(:,:)
   enddo  !  ie=nets,nete
#ifdef OPENACC_HOMME
!$acc end parallel loop 
!!$acc update device(elem,hvcoord)
!$acc parallel loop gang vector collapse(2) present(elem,dp_ie,dp_star_ie,hvcoord) 
#endif
   do ie=nets,nete         
     do k=1,nlev
        dp_ie(ie,:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
        if (rsplit==0) then
           dp_star_ie(ie,:,:,k) = dp_ie(ie,:,:,k) + dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) -&
                elem(ie)%derived%eta_dot_dpdn(:,:,k))
        else
           dp_star_ie(ie,:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)
        endif
     enddo
   enddo  !  ie=nets,nete
#ifdef OPENACC_HOMME
!$acc end parallel loop
!$acc update self(dp_star_ie)
#endif  
   do ie=nets,nete   
     if (minval(dp_star_ie(ie,:,:,:))<0) then
#ifdef OPENACC_HOMME
!$acc update self(dp_ie,elem)
#endif     
        do k=1,nlev
        do i=1,np
        do j=1,np
           if (dp_star_ie(ie,i,j,k ) < 0) then
              print *,'index ie,i,j,k = ',ie,i,j,k
              print *,'dp_star = ',dp_star_ie(ie,i,j,k)
              !print *,'dp,dp_star = ',dp_ie(ie,i,j,k),dp_star_ie(ie,i,j,k)
              !print *,'eta_dot_dpdn = ',elem(ie)%derived%eta_dot_dpdn(i,j,k+1),elem(ie)%derived%eta_dot_dpdn(i,j,k)
              !print *,"column location lat,lon (radians):",elem(ie)%spherep(i,j)%lat,elem(ie)%spherep(i,j)%lon
           endif
        enddo
        enddo
        enddo
        call abortmp('negative layer thickness.  timestep or remap time too large')
     endif
   enddo  !  ie=nets,nete

     if (rsplit>0) then
        !removing theta_ref does not help much and will not conserve theta*dp
        !call set_theta_ref(hvcoord,dp_star,theta_ref)

       call get_phinh_openacc(hvcoord,elem,np1,dp_star_ie,phi_ref_ie,nets,nete) 

#ifdef OPENACC_HOMME
!$acc parallel loop gang present(elem,dp_star_ie,phi_ref_ie,ttmp_ie)
#endif    
       do ie=nets,nete
         elem(ie)%state%phinh_i(:,:,:,np1)=&
             elem(ie)%state%phinh_i(:,:,:,np1) -phi_ref_ie(ie,:,:,:)
         !  REMAP u,v,T from levels in dp3d() to REF levels
         ttmp_ie(ie,:,:,:,1)=elem(ie)%state%v(:,:,1,:,np1)*dp_star_ie(ie,:,:,:)
         ttmp_ie(ie,:,:,:,2)=elem(ie)%state%v(:,:,2,:,np1)*dp_star_ie(ie,:,:,:)
         ttmp_ie(ie,:,:,:,3)=elem(ie)%state%vtheta_dp(:,:,:,np1)   ! - theta_ref*dp_star*Cp
         do k=1,nlev
            ttmp_ie(ie,:,:,k,4)=elem(ie)%state%phinh_i(:,:,k+1,np1)-&
                elem(ie)%state%phinh_i(:,:,k,np1) 
            ttmp_ie(ie,:,:,k,5)=elem(ie)%state%w_i(:,:,k+1,np1)-&
                elem(ie)%state%w_i(:,:,k,np1)
         enddo        
       enddo  !  ie=nets,nete
#ifdef OPENACC_HOMME
!$acc end parallel loop
!$acc update self(ttmp_ie,dp_star_ie, dp_ie)
#endif         

       call remap1_1_openacc(ttmp_ie,np,5,dp_star_ie,dp_ie,nets,nete)
       
#ifdef OPENACC_HOMME
!$acc parallel loop gang present(elem,ttmp_ie, dp_ie)
#endif        
       do ie=nets,nete 
         !call set_theta_ref(hvcoord,dp,theta_ref)
         elem(ie)%state%v(:,:,1,:,np1)=ttmp_ie(ie,:,:,:,1)/dp_ie(ie,:,:,:)
         elem(ie)%state%v(:,:,2,:,np1)=ttmp_ie(ie,:,:,:,2)/dp_ie(ie,:,:,:)
         elem(ie)%state%vtheta_dp(:,:,:,np1)=ttmp_ie(ie,:,:,:,3) ! + theta_ref*dp*Cp
                 
         do k=nlev,1,-1
            elem(ie)%state%phinh_i(:,:,k,np1)=&
                elem(ie)%state%phinh_i(:,:,k+1,np1)-ttmp_ie(ie,:,:,k,4)  !/dp(:,:,k)
            elem(ie)%state%w_i(:,:,k,np1)=&
                elem(ie)%state%w_i(:,:,k+1,np1)-ttmp_ie(ie,:,:,k,5)  !/dp(:,:,k)
         enddo
       enddo  !  ie=nets,nete
#ifdef OPENACC_HOMME
!$acc end parallel loop
#endif  
       
       call get_phinh_openacc(hvcoord,elem,np1,dp_ie,phi_ref_ie,nets,nete)

#ifdef OPENACC_HOMME
!$acc parallel loop gang present(elem,phi_ref_ie)
#endif        
       do ie=nets,nete  
         elem(ie)%state%phinh_i(:,:,:,np1)=&
             elem(ie)%state%phinh_i(:,:,:,np1)+phi_ref_ie(ie,:,:,:)
             
         ! since u changed, update w b.c.:
         elem(ie)%state%w_i(:,:,nlevp,np1) = (elem(ie)%state%v(:,:,1,nlev,np1)*elem(ie)%derived%gradphis(:,:,1) + &
             elem(ie)%state%v(:,:,2,nlev,np1)*elem(ie)%derived%gradphis(:,:,2))/g
       
       enddo  !  ie=nets,nete
#ifdef OPENACC_HOMME
!$acc end parallel loop
#endif        
     endif
     ! remap the gll tracers from lagrangian levels (dp_star)  to REF levels dp
     if (qsize>0) then
       
        call remap1_3_openacc(elem,np1_qdp,np,qsize,dp_star_ie,dp_ie,nets,nete)

     endif
#ifdef OPENACC_HOMME
!!$acc update self(elem)
!$acc exit data delete(dp_ie,dp_star_ie, phi_ref_ie, ttmp_ie)
#endif   
  call t_stopf('vertical_remap')
  end subroutine vertical_remap


end module 




