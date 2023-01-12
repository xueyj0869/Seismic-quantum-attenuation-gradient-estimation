% Sample code of the paper:
% 
% Ya¨Cjuan Xue, Xing-jian Wang, Jun-xing Cao, Zhe-ge Liu and Jia Yang,
% "Quantum mechanics¨Cbased seismic energy absorption analysis for hydrocarbon detection." 
% submitted to Geophysical Journal International (2023).
% 
% MATLAB code prepard by Ya-juan Xue
% E-mail: xueyj0869@163.com
% 
% This function is our quantum attenuation gradient estimation algorithm 

% During this function, we call the subfunction of the Hamiltonian associates with the signal
% and the associated eigenvalues for signal using Quantum adaptative basis(QAB)
% prepared by Sayantan Dutta et al."Quantum mechanics-based signal and image denoising." 
% arXiv preprint arXiv:2004.01078 (2020).

function out=quanattenver(seis,fs,h)
warning('off')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% seis:input seismic data;
% fs:sampling frequency
% h: Planck constant;   eg. h=0.016;

% Output
% out: quantum attenuation gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[trace,time]=size(seis);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:trace
    i
   s=seis(i,:);
   S=s./max(s);
   [psi,E] = f_ondes1D(S,h);
   E=E';   
   for j=1:size(psi,2)
%        j
       ss=psi(:,j);
       y=fft(ss);
       mag=abs(y.^2);
       smag=mag(1:round(time/2));
       slga=(log(smag+1));
       fq=sqrt(2*(E-S(j)))./h;   
       fqr=real(fq);
      ind=find(slga==max(slga));
      ind=ind(1);  
      ind2=ind+10;
      if ind2>time/2
          ind22=time/2;
      else
          ind22=ind2;
      end
      fx=fq(ind:round(ind22));
      fab=slga(ind:round(ind22))';

       kk=polyfit(fx,fab,1);
       afirst=kk(1);

       out(i,j)=-real(afirst);
   end 
end

