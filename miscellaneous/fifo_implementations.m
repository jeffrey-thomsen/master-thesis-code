%% FIFO implementations
% https://de.mathworks.com/matlabcentral/answers/289883-fast-fifo-array-other-datatype

%% Walter Roberson 2016
array(1:end-1) = array(2:end);
array{end} = newvalue;

%% Andrew Ward 2018
Nmax=50000;
A1=1:Nmax;
Ap=A1;


%% 
%Method 1 concatenate array with element - Slow
tic
for i=1:Nmax
    
    A1=[A1(2:end),Nmax+i];
    
end
toc
%% 


%% Method 2 overwite ellemnts followed by a circular shift - Fastest - 
% avoids taking a sub portion of the array
% which I think is costly, just replaces the value, 1,2,3,4 
% becomes 5,2,3,4, so to get it in order you need to
% do the circular shift at the end.
tic

A2=1:Nmax;

for i=1:Nmax
    ic=mod(i,Nmax+1);
    A2(ic)=Nmax+i;
end
A2=circshift(A2,-i);
toc

sum(A2~=A1)

%% Method from above, still slow
A3=Ap;
tic
for i=1:Nmax
    
    A3(1:end-1)=A3(2:end);
    A3(end)=Nmax+i;
end

toc
sum(A3~=A1);