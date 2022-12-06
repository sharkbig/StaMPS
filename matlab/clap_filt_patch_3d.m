function [phout]=clap_filt_patch(ph,alpha,beta,low_pass)
%CLAP_FILT_PATCH Combined Low-pass Adaptive Phase filtering on 1 patch
%   [ph_out]=clap_filt_patch(ph,alpha,beta)
%
%
%   Andy Hooper, June 2006
%  codegen clap_filt_patch_3d -args {coder.typeof(single(1j),[inf inf inf]),0.5,0.1,coder.typeof(0,[inf])}
B=gausswin(7)*gausswin(7)';
ph(isnan(ph))=0;

ph_fft=fft(ph,[],1);
ph_fft=fft(ph_fft,[],2);

H=abs(ph_fft);
H=fftshift(H,1);
H=fftshift(H,2); 

% H=filter2_3d_mex(B,H);
H=convn(H,B,'same');

H=ifftshift(H,1);
H=ifftshift(H,2);

meanH=squeeze(median(H,[1 2]));
meanH(meanH==0)=1;
[a b c]=size(H);
meanH(meanH==0)=1;
meanH=repmat(meanH,1,1,a*b);
meanH=reshape(permute(meanH,[3 1 2]),[a b c]);
H=H./meanH;

H=H.^alpha;
H=H-1;
H(H<0)=0;
G=H.*beta+low_pass;

phout=ph_fft.*G;   
phout=ifft(phout,[],1);
phout=ifft(phout,[],2); 
