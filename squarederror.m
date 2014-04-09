function sse=squarederror(params,funchandle,Input,Actual_Output)
%sum squared error evaluation function for fminsearch routine
%using anonymous function
%
%Saul Kato

Fitted_Curve=funchandle(params,Input);
Error_Vector=Fitted_Curve - Actual_Output;
sse=sum(Error_Vector.^2);

end