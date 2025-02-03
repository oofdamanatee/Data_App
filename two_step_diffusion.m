function out = two_step_diffusion(t,D1,D2,A,C1,C2,nmax,sigma,dx1,dx2,dy,rlim)
    out = diffusion_moving_beam(t,D1,A,C1,nmax,sigma,dx1,dy,'rlim',rlim)...
        + diffusion_moving_beam(t,D2,A,C2,nmax,sigma,dx2,dy,'rlim',rlim);
end