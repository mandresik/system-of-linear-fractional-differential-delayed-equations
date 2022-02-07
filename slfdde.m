function [t,x]=slfdde(met,limit,t0,T,tau,h,alpha,x0,A,B,f)
%% INPUTS DESCRIPTION
%   met   ... method used for numerical integration
%              met=0 for explicit rectangular method 
%              met=1 for implicit rectangular method
%   t0    ... initial time
%   T     ... ending time
%   tau   ... delay
%   h     ... step
%   alpha ... order
%   x0    ... initial function 
%   A     ... matrix A         
%   B     ... matrix B          
%   f     ... vector f  
%
%   Inputs example
%   met=0;limit=100;t0=0;T=70;tau=pi;h=0.01;alpha=2/5;
%   x0=@(t)[sin(t)*cos(t);sin(t)*cos(t);(cos(t))^2-(sin(t))^2;(cos(t))^2-(sin(t))^2];
%   A=@(t)[0 0 1 0;0 0 0 1;0 -2 0 0;-2 0 0 0];
%   B=@(t)[0 0 0 0;0 0 0 0;-2 0 0 0;0 -2 0 0];
%   f=@(t)[0;0;0;0]; 
%% INPUTS VERIFICATION
m=ceil(alpha); % ceil func of alpha used for computation
var_num=length(A(0)); % number of variables
if (m~=size(x0(0),2)) % verification of number of initial conditions 
    if (m>size(x0(0),2))
        msgbox(sprintf('                         Not enough initial conditions assigned,\nnumber of initial conditions must be first integer greater or equal than alpha.\n\n                       alpha = %.2f ==> number of conditions = %d \n',alpha,m),'Error - conditions')
    else
        msgbox(sprintf('                          Too much initial conditions assigned,\nnumber of initial conditions must be first integer greater or equal than alpha.\n\n                       alfa = %.2f ==> number of conditions = %d \n',alpha,m),'Error - conditions')
    end
    return
end % verification of dimension
if(size(A(0),1)~=var_num || size(A(0),2)~=var_num || size(B(0),1)~=var_num || size(B(0),2)~=var_num || size(x0(0),1)~=var_num || size(f(0),1)~=var_num)
    msgbox(sprintf('Please check dimensions of matrices and vectors,\n dimensions of matrices A and B must be %d x %d,\n dimensions of vectors x0 and f must be %d x 1.',var_num,var_num,var_num),'Error - dimensions')
    return 
end % message about default method choice
if(met~=0 && met~=1)
    msgbox({'Number met for method choice is neither 0 and 1,', 'for computation is used explicit method.'},'Method choice')
end
%% COMPUTATION
t_tau=t0-tau; % time serie t creation
temp=t0-h:-h:t_tau;
t=[temp(end:-1:1) t0:h:T];

index=length(temp);
N=length(t)-length(temp)-1;

x=zeros(var_num,length(t)); % alocation of x and b
b=zeros(1,N); 
F=zeros(var_num,N); % store of matrices multiplication
x0TP=x0(t0); % store of init cond for extracting cols in Taylor polynomial

for i=1:index+1 % storing init conditions to x in time leq t0
    x(:,i)=x0(t(i))*[1;zeros(m-1,1)]; % first col in x0
end

At=function_values(A,limit);
Bt=function_values(B,limit);
ft=function_values(f,limit);

% method call
if(met==1)
    implicit
else
    explicit
end

function Mt=function_values(M,limit)
    cols=size(M(0),2); % number of columns
    risk_index=zeros(var_num,1); % vector of risk indexes
    col=0; % index for saving risk indexes
    M0=M(0);
    for r=1:var_num % indexes that are +-infinity fot t=0
        for s=1:cols
            if(M0(r,s)==Inf || M0(r,s)==-Inf)
            risk_index=[risk_index(:,1:col) [r;s]];
            col=col+1;
            end
        end
    end

    Mt=zeros(var_num,cols*length(t));

    for j=1:length(t) % function values
        Mt(:,cols*j-cols+1:cols*j)=M(t(j));
    end

    if(risk_index~=zeros(var_num,1)) % if there is risk_index, vaules of t between -1 and 1 are investigated (possibly set to +-limit) 
    risk_start_index=sum(t<=-1); % index that multiplicates risk_index
    inv_num=sum(t<1)-sum(t<=-1);
    risks_num=size(risk_index,2);
    for j=1:inv_num
        for k=1:risks_num
            if(Mt(risk_index(1,k),risk_start_index*cols+risk_index(2,k)+cols*(j-1))>limit)
                Mt(risk_index(1,k),risk_start_index*cols+risk_index(2,k)+cols*(j-1))=limit;
            elseif (Mt(risk_index(1,k),risk_start_index*cols+risk_index(2,k)+cols*(j-1))<-limit)
                Mt(risk_index(1,k),risk_start_index*cols+risk_index(2,k)+cols*(j-1))=-limit;
            end
        end
    end
    end
end
    
% explicit method algorithm
function explicit
    for n=1:N 
        col_num=size(A(0),2); 
        ind=col_num*(n+index); 
        b(n)=(n^alpha-(n-1)^alpha)/gamma(alpha+1);
        F(:,n)=At(:,ind-col_num+1:ind)*x(:,n+index)+Bt(:,ind-col_num+1:ind)*x(:,n)+ft(:,n+index);
        TP=zeros(var_num,1);
        for k=0:m-1
            TP=TP+((t(n+index+1)-t0)^k).*x0TP(:,k+1)./factorial(k);
        end                   
        sum=zeros(var_num,1);
        for j=1:n
            sum=sum+b(n-j+1).*F(:,j);
        end
        x(:,n+index+1)=TP+h^alpha.*sum;
    end
end

% implicit method algorithm
function implicit
    sum=zeros(var_num,1);
    for n=1:N
        col_num=size(A(0),2); 
        ind=col_num*(n+index+1); 
        TP=zeros(var_num,1);
        for k=0:m-1
            TP=TP+((t(n+index+1)-t0)^k).*x0TP(:,k+1)./factorial(k); 
        end                   
        x(:,n+index+1)=(eye(var_num)-h^alpha/gamma(alpha+1)*At(:,ind-col_num+1:ind))\(TP+h^alpha/gamma(alpha+1)*(Bt(:,ind-col_num+1:ind)*x(:,n+1)+f(t(n+index+1))+sum)); 
        b(n)=(n+1)^alpha-n^alpha;
        F(:,n)=At(:,ind-col_num+1:ind)*x(:,n+index+1)+Bt(:,ind-col_num+1:ind)*x(:,n+1)+ft(n+index+1);
    
        sum=zeros(var_num,1);
        for j=1:n 
            sum=sum+b(n-j+1).*F(:,j);
        end 
    end
end

end