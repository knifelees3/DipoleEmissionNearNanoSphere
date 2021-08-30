function nP=Fun_nP_DipoleNearSphere(WL0,nb,ns,RAu,dis)
	% Only the total decay rate will be simulated
	k0=2*pi/WL0;
	kb=k0*nb;
	ks=k0*ns;
	kbRAu=kb*RAu;
	ksRAu=ks*RAu;
	[num_dis,~]=size(dis(:));
	num_sum=50; % default number of modes to be sumed
	nP=zeros(num_dis,2); % with nP(:,1) store the np_ortho and nP(:,2) store the np_para
	for m=1:num_dis
		rate_ortho=0;
		rate_para=0;
		r=dis(m);
		for l=1:num_sum
			a1=(ns^2*spBJ(l,ksRAu)*dspBJ(l,kbRAu)-spBJ(l,kbRAu)*dspBJ(l,ksRAu));
			a2=(ns^2*spBJ(l,ksRAu)*dspBH(l,kbRAu)-spBH(l,kbRAu)*dspBJ(l,ksRAu));   
			b1=(spBJ(l,ksRAu)*dspBJ(l,kbRAu)-spBJ(l,kbRAu)*dspBJ(l,ksRAu));
			b2=(spBJ(l,ksRAu)*dspBH(l,kbRAu)-spBH(l,kbRAu)*dspBJ(l,ksRAu));
			x=kb*r;
			al=a1/a2;
			bl=b1/b2;
			rate_ortho=rate_ortho+(2*l+1)*l*(l+1)*(-al)*(spBH(l,x)/(x))^2;
			rate_para=rate_para+(l+0.5)*((-al)*(dspBH(l,x)/(x))^2+(-bl)*(spBH(l,x))^2);   
		end
		rate_ortho=real(rate_ortho)*3/2+1;
		rate_para=1+1.5*real(rate_para);
		nP(m,1)=rate_ortho;
		nP(m,2)=rate_para; 
	end
end