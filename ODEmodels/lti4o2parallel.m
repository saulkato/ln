function outp=lti4o2parallel(pm,t)


% odeModel.index=3;
% odeModel.name='4o2parallel';
% odeModel.numParams=6;
% odeModel.paramLabels={'TriseF','TriseS','TfallF','TfallS','kas/kaf','kas+kaf'};
% odeModel.func=@lti4o2parallelgc;
% odeModel.funcDeconvolved=@lti4o2parallel;
% odeModel.paramMax=[30,30,30,40,100,1];
% odeModel.paramMin=[0.001,0.001,0.001,0.001,.00001,-1];
% odeModel.standardLossFuncFlag=true;
% odeModel.fitseed=[20 10 1.051 37 .026 -.019];

tau_s_rise=pm(1);
tau_f_rise=pm(2);
tau_s_fall=pm(3);
tau_f_fall=pm(4);

kasoverkaf=pm(5);
kaspluskaf=pm(6);

kaf=kaspluskaf/(kasoverkaf+1);

kas=kaspluskaf*(kasoverkaf)/(kasoverkaf+1);


outp=kaf*(tau_f_rise*tau_f_fall)*(1-exp(-t/tau_f_rise)).*exp(-t/tau_f_fall)/abs(tau_f_rise-tau_f_fall) + ...
     kas*(tau_s_rise*tau_s_fall)*(1-exp(-t/tau_s_rise)).*exp(-t/tau_s_fall)/abs(tau_s_rise-tau_s_fall);
 
end