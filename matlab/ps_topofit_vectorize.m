function [K0,C0,coh0,phase_residual]=ps_topofit_vectorize(cpxphase,bperp,n_trial_wraps,asym)
%PS_TOPOFIT find best-fitting range error 
%   PS_TOPOFIT(cpxphase,bperp,n_trial_wraps,plotflag,asym)
%   ASYM = -1 (only -ve K searched) to +1 (only +ve K) (default 0)
%
%   Andy Hooper, June 2006
%
%   ==========================================================
%   04/2007 AH: Added 64-bit machine compatibility
%   04/2007 AH: Tightened up max topo error processing
%   11/2012 AH: Add asymmetry option
%   ==========================================================
% 
%  Jun-Yan Chen
%   08/2022 JY: make the code vectorize to speed up,
%               remove plot flag
    

if nargin<4
    asym=0;
end

n_ifg=size(cpxphase,2);
bk_size=size(cpxphase,1);

bperp_range=max(bperp,[],2)-min(bperp,[],2);
trial_mult=[-ceil(8*n_trial_wraps):ceil(8*n_trial_wraps)]+asym*8*n_trial_wraps;
n_trials=length(trial_mult);
trial_phase=pi/4*bperp.*(1./bperp_range);
cpxmat=exp(kron(-j*trial_phase,trial_mult));
cpxmat=reshape(transpose(cpxmat),n_trials,n_ifg,bk_size);
cpxmat=permute(cpxmat,[3 2 1]);
phaser=cpxmat.*cpxphase;
phaser_sum=squeeze(sum(phaser,2));

coh_trial=abs(phaser_sum)./sum(abs(cpxphase),2);
[dummy,coh_max_ix]=max(coh_trial,[],2);

K0=pi/4./bperp_range.*trial_mult(coh_max_ix)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% linearise and solve %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

resphase=cpxphase.*exp(-j*(K0.*bperp));
offset_phase=sum(resphase,2);
resphase=angle(resphase.*conj(offset_phase));
weighting=abs(cpxphase);
%%%% no piese-wise inversion function %%%%
for k=1:bk_size
    mopt=double(weighting(k,:).*bperp(k,:))'\double(weighting(k,:).*resphase(k,:))';
    K0(k)=K0(k)+mopt;
end

phase_residual=cpxphase.*exp(-j*(K0.*bperp));
mean_phase_residual=sum(phase_residual,2);
C0=angle(mean_phase_residual);
coh0=abs(mean_phase_residual)./sum(abs(phase_residual),2);