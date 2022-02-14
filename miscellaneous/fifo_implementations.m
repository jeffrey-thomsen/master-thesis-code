%% My reference experiments
% All equal results, Method 1 more than twice as fast as Methods 2 and 3

% testSS = complex Nx1 signal array
% testFIFO = 1xM array

% Method 1: use circshift
tic
testFIFO1=testFIFO;
for i=1:length(testSS)
    testFIFO1(end) = testSS(i);
    testFIFO1 = circshift(testFIFO1,1);
end
toc

% Method 2: concatenate
tic
testFIFO2=testFIFO;
for i=1:length(testSS)
    testFIFO2 = [testSS(i), testFIFO2(1:end-1)];
end
toc

% Method 3: just overwrite entire array
tic
testFIFO3=testFIFO;
for i=1:length(testSS)
    testFIFO3(2:end) = testFIFO3(1:end-1);
    testFIFO3(1) = testSS(i);
end
toc

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