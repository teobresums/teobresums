function df=Deriv4(f,x,order)
%-----------------------------------------------
% Computes the first derivative of the function
% Centered but at the edges.
% USAGE: df = Deriv4(f,x,order) where order can
%             be 1,2,4
%-----------------------------------------------

%disp('>> INFO: EOB_D1>>');
%disp(['Computing first derivative at ',num2str(order),' order']);
Nmin = 1;   
Nmax = length(f);
% preallocation for increasing speed!
% Note that this is 1 columns and N rows
df = zeros(Nmax,1);

if order==2
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %% 2nd order centered stancil
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=Nmin+1:Nmax-1
        df(i) = (f(i+1)-f(i-1))/(x(i+1)-x(i-1));
    end 
    
    %% 2nd order boundaries
    df(Nmin) = -(3*f(Nmin)-4*f(Nmin+1)+f(Nmin+2))/(x(Nmin+2)-x(Nmin));  %+o(dx^2)
    df(Nmax) =  (3*f(Nmax)-4*f(Nmax-1)+f(Nmax-2))/(x(Nmax)-x(Nmax-2));  %+o(dx^2)
    
    
elseif order==4
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %% 4th order centered stancil
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    %for i=Nmin+2:Nmax-2
    %   df(i) = 1/3*( 8*f(i+1)-f(i+2) - 8*f(i-1)+f(i-2))/(x(i+2)-x(i-2));
    %end

    fim2 = f(1:end-4);
    fim1 = f(2:end-3);
    fip1 = f(4:end-1);
    fip2 = f(5:end);
    xip2 = x(5:end);
    xim2 = x(1:end-4);
    dx   = xip2 - xim2;
    
    df(Nmin+2:Nmax-2) = 1/3*(8*fip1 - fip2 - 8*fim1 + fim2)./dx;               
    
    % 4th order boundaries corrected from Mathematica sheet
    df(Nmin)   = (-25*f(Nmin) + 48*f(Nmin+1) - 36*f(Nmin+2) +16*f(Nmin+3)-3*f(Nmin+4))/(3*(x(Nmin+4)-x(Nmin)));
    df(Nmin+1) = (-3*f(Nmin)-10*f(Nmin+1)+18*f(Nmin+2)-6*f(Nmin+3)+f(Nmin+4))/(3*(x(Nmin+4)-x(Nmin)));
    
    df(Nmax)   = (25*f(Nmax) - 48*f(Nmax-1) + 36*f(Nmax-2) - 16*f(Nmax-3)+3*f(Nmax-4))/(3*(x(Nmax)-x(Nmax-4)));
    df(Nmax-1) = (-f(Nmax-4)+6*f(Nmax-3)-18*f(Nmax-2)+10*f(Nmax-1)+3*f(Nmax))/(3*(x(Nmax)-x(Nmax-4)));        
  
else
    error('******************* "order" must be 2 or 4 ************************')      
end

return