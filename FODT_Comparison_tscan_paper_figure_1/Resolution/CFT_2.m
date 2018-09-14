% Modified to start at 0 frequency to get phase starting at 0
% Remove periodic component by only considering points below F(M)

function [w,F]=CFT_2(t,f)

N=length(t);			% # of sampled points
dt=t(2)-t(1);			% temporal resolution; temporal window = N*dt, not (N-1)*dt
dw=2*pi/N/dt;           % spectral resolution

F=fft(f);               % F[k]=DFT{f[n]}; f[n]=f(t=t0+n*dt), k,n=0,1,...,(N-1)
                        % fftshift([1 2, 3])=[3, 1 2]; w=0 at center
                        % fftshift([1 2, 3 4])=[3 4, 1 2]; w=0 is biased right

if mod(N,2)==0			% if N is even
   M=N/2;               % G[M]=F[0]~F(w=0), if G[q]=fftshift{F[k]}, q=0,1,...N-1
   w=[0:2*(N-M)-1]*dw;      % sampled angular frequencies: p*dw
else					% if N is odd
   M=(N-1)/2;
   w=[0:2*(N-M-1)]*dw;      % sampled angular frequencies: p*dw
end

                        
%F(M+1:2*(N-M-1)+1)=0; 
F=F*dt;
%F=fftshift(F)*dt.*exp(-j*w*t(1));