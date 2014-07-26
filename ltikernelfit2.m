function [outp,outptype,outpm]=ltikernelfit2(kern,modelNameOrFunc,dt,starting_paramguess,lnmstruct,standardLossFuncFlag)
%[outp,outptype,outpm]=ltikernelfit2(kern,type,dt,starting_paramguess)
%fit analytical kernels to measured kernel, via lneval
%update 20120917 variable starting params
%Saul Kato
%
if (nargin<6) standardLossFuncFlag=true; end
if (nargin<3) dt=0.05; end

inp=dt*(0:(length(kern)-1))';

options=optimset('Display','off');  %'iter' 'final' 'none'

      
if standardLossFuncFlag
        
        [pm,fval]=fminsearch(@squarederror,starting_paramguess,options,modelNameOrFunc,inp,kern);
        outp=modelNameOrFunc(pm,inp);
        outptype=modelNameOrFunc;
        outpm=pm;
      
else
    
        if strcmp(modelNameOrFunc,'over')
            %Starting=[0.0202    0.8066   0.0030];
            if nargin<4
                starting_paramguess=[    0.1714    3.9275   0.0477];
            end
            [pm_over,fval]=fminsearch(@overddifffit,starting_paramguess,options,inp,kern);
            outp=overddifffn(pm_over,inp);
            outptype=modelNameOrFunc;
            outpm=pm_over;

        elseif strcmp(modelNameOrFunc,'over2')
            if nargin<4
                starting_paramguess=[-.5,.3,-.6];
            end
            [pm,fval]=fminsearch(@lti2ofit,starting_paramguess,options,inp,kern);
            outp=lti2ofn(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;     

        elseif strcmp(modelNameOrFunc,'under')
            if nargin<4
                starting_paramguess=[1.28   0.24    0.0426];
            end
            %[pm_under,fval]=fminsearch(@underdampedfit,Starting,options,inp,kern);
            %outp=underdampedfn(pm_under,inp);
            %z=Starting; 
            %g=.3*sqrt( ((z(2)/z(3))^2)/(1+(z(2)/z(3))^2) );    
            %pnew=[ 2*g/z(2)    .01*g     -z(1)/z(3)] 
            %Starting=pnew;
            [pm_under,fval]=fminsearch(@underddifffit,starting_paramguess,options,inp,kern);
            %pm_under=Starting;
            outp=underddifffn(pm_under,inp);

            outptype=modelNameOrFunc;
            outpm=pm_under;

        elseif strcmp(modelNameOrFunc,'crit')
            if nargin<4
                starting_paramguess=[1.9675   -0.1038];
            end
            [pm,fval]=fminsearch(@alphadifffit,starting_paramguess,options,inp,kern);
            outp=alphadifffn(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm; 

        elseif strcmp(modelNameOrFunc,'3offc1')  %feed-forward cascade type one (non-branching)
            if nargin<4
                starting_paramguess=[-1 -2 -3 5];
            end
            [pm,fval]=fminsearch(@lti3offc1fit,starting_paramguess,options,inp,kern);
            outp=lti3offc1fn(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm; 

        elseif strcmp(modelNameOrFunc,'3o1')   
            if nargin<4
                starting_paramguess=[-.5,-.1,.3,0,-0.3,-2];
            end
            [pm,fval]=fminsearch(@lti3ofit,starting_paramguess,options,inp,kern);
            outp=lti3ofn(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;     

        elseif strcmp(modelNameOrFunc,'3o3')
            if nargin<4
                starting_paramguess=[-.5,-.1,.3,0,-0.3];
            end
            [pm,fval]=fminsearch(@lti3ofit3,starting_paramguess,options,inp,kern);
            outp=lti3ofn3(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;     

        elseif strcmp(modelNameOrFunc,'3o4')
            if nargin<4
                 starting_paramguess=[2 2 20 -5 -15];
            end
            [pm,fval]=fminsearch(@lti3ofit4,starting_paramguess,options,inp,kern);
            outp=lti3ofn4(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;

        elseif strcmp(modelNameOrFunc,'3o4gc')
            if nargin<4
                 starting_paramguess=[2 2 20 -5 -15];
            end
            [pm,fval]=fminsearch(@lti3odifffit4,starting_paramguess,options,inp,kern);
            outp=lti3odifffn4(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;

        elseif strcmp(modelNameOrFunc,'3o5gc')
            if nargin<4
                 starting_paramguess=[2 2 20 -5 -15];
            end
            [pm,fval]=fminsearch(@lti3odifffit5,starting_paramguess,options,inp,kern);
            outp=lti3odifffn5(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;

        elseif strcmp(modelNameOrFunc,'3o5wang')
            if nargin<4
                 starting_paramguess=[2 2 20 -5 -15];
            end
            [pm,fval]=fminsearch(@lti3odifffitwang,starting_paramguess,options,inp,kern);
            outp=lti3odifffnwang(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;

        elseif strcmp(modelNameOrFunc,'3o5wang2')
            if nargin<4
                 starting_paramguess=[2 2 20 -5 -15];
            end
            [pm,fval]=fminsearch(@lti3odifffitwang2,starting_paramguess,options,inp,kern);
            outp=lti3odifffnwang(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;

        elseif strcmp(modelNameOrFunc,'3o5wangp2')
            disp(['fitting ' modelNameOrFunc]);
            if nargin<4
                 starting_paramguess=[2 2 20 -5 -15];
            end
            [pm,fval]=fminsearch(@lti3odifffitwangp2,starting_paramguess,options,inp,kern);
            outp=lti3odifffnwangp2(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;

        elseif strcmp(modelNameOrFunc,'3o5wangp2_lneval2')
            disp(['fitting ' modelNameOrFunc]);
            if nargin<4
                 starting_paramguess=[2 2 20 -5 -15];
            end
            [pm,fval]=fminsearch(@lti3odifffitwangp2_lneval2,starting_paramguess,options,lnmstruct);
            outp=lti3odifffnwangp2(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;

        elseif strcmp(modelNameOrFunc,'3o5gcb')
            if nargin<4
                 starting_paramguess=[2 2 20 -5 -15];
            end
            [pm,fval]=fminsearch(@lti3odifffit5b,starting_paramguess,options,inp,kern);
            outp=lti3odifffn5b(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;


        elseif strcmp(modelNameOrFunc,'3o5gcexact')
            if nargin<4
                 starting_paramguess=[2 2 20 -5 -15];
                 disp('using default 5-value starting paramguess');
            end
            %size(starting_paramguess)
            [pm,fval]=fminsearch(@lti3odifffit5exact,starting_paramguess,options,inp,kern);
            outp=lti3odifffn5exact(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;      

        elseif strcmp(modelNameOrFunc,'3o5gcexactneg')
            if nargin<4
                 starting_paramguess=[2 2 20 -5 -15];
                 disp('using default 5-value starting paramguess');
            end
            size(starting_paramguess)
            [pm,fval]=fminsearch(@lti3odifffit5exactneg,starting_paramguess,options,inp,kern);
            outp=lti3odifffn5exactneg(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;      

        elseif strcmp(modelNameOrFunc,'3oo2')  %Manuel's O2 chip
            if nargin<4
                 starting_paramguess=[2 2 20 -5 -15];
            end
            [pm,fval]=fminsearch(@lti3oo2fit,starting_paramguess,options,inp,kern);
            outp=lti3oo2fn(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;

        elseif strcmp(modelNameOrFunc,'3oo22')  %Manuel's O2 chip V2
            if nargin<4
                 starting_paramguess=[2 2 20 -5 -15];
            end
            [pm,fval]=fminsearch(@lti3oo2fit2,starting_paramguess,options,inp,kern);
            outp=lti3oo2fn2(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;

        elseif strcmp(modelNameOrFunc,'4o1')
            if nargin<4
                 starting_paramguess=[2 2 2 20 -5 -15];
            end
            [pm,fval]=fminsearch(@lti4ofit1,starting_paramguess,options,inp,kern);
            outp=lti4ofn1(pm,inp);
            outptype=modelNameOrFunc;
            outpm=pm;

        elseif strcmp(modelNameOrFunc,'best')

            starting_paramguess=[    0.1714    3.9275   0.0477];
            [pm_over,fval]=fminsearch(@overddifffit,starting_paramguess,options,inp,kern);
            outp_o=overddifffn(pm_over,inp);
            disp(['overd: ' num2str(vaf(outp_o,kern))]);

            starting_paramguess=[-0.0050    0.1294    0.3226];
            [pm_under,fval]=fminsearch(@underddifffit,starting_paramguess,options,inp,kern);
            outp_u=underddifffn(pm_under,inp);
            disp(['underd: ' num2str(vaf(outp_u,kern))]);

            starting_paramguess=[1.9675   -0.1038];
            [pm,fval]=fminsearch(@alphadifffit,starting_paramguess,options,inp,kern);
            outp_a=alphadifffn(inp,pm(2),pm(1));
            disp(['crit: ' num2str(vaf(outp_a,kern))]);

            if vaf(outp_o,kern) > vaf(outp_u,kern)
                if vaf(outp_o,kern)>vaf(outp_a,kern)
                    outp=outp_o;
                    outptype='over';
                    outpm=pm_over;
                else
                    outp=outp_a;
                    outptype='crit';   
                    outpm=pm;
                end
            elseif vaf(outp_u,kern) > vaf(outp_a,kern)
                    outp=outp_u;
                    outptype='under';
                    outpm=pm_under;
            else
                    outp=outp_a;
                    outptype='crit';
                    outpm=pm;
            end


        end

end