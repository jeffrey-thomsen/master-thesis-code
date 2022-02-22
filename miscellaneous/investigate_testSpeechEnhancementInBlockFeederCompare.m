% run this for investigation of testSpeechEnhancementInBlockFeederCompare
% during a breakpoint before verification

% compare output signals
numel(find(abs(batchOutput-blockOutput)<1e-16))

% compare enhanced subband signals
numel(find(abs(batchESB.L-blockESB.L)>1e-9))
[rows,cols]=(find(abs(batchESB.L-blockESB.L)>1e-9));
% see which subband samples are not identical
scatter(rows,cols)

% compare IVS mask logical arrays
for i=1:500;for j=1:30;batchIvsMatrix(i,j)=batchivsMask{j}(i);end;end;
for i=1:500;for j=1:30;blockIvsMatrix(i,j)=blockivsMask{i}{j};end;end;
numel(find(batchIvsMatrix~=blockIvsMatrix))

% compare computed IPD values
for i=1:500;for j=1:30;batchIpdRadMatrix(i,j)=batchipdRad{j}(i);end;end;
for i=1:500;for j=1:30;blockIpdRadMatrix(i,j)=blockipdRad{i}{j};end;end;
numel(find(abs(batchIpdRadMatrix-blockIpdRadMatrix)<1e-14))

% compare number of samples with passed azimuth (integer)
blockAzSum=0;for i=1:500;for j=1:30;blockAzSum=blockAzSum+numel(blockazimuthRadCells{i}{j});end;end;
batchAzSum=0;for i=1:30;batchAzSum=batchAzSum+numel(batchazimuthRadCells{i});end;

% compare cumulative sum of passed azimuth values (numeric/float)
blockAzSumVal=0;for i=1:500;for j=1:30;blockAzSumVal=blockAzSumVal+sum(blockazimuthRadCells{i}{j});end;end;
batchAzSumVal=0;for i=1:30;batchAzSumVal=batchAzSumVal+sum(batchazimuthRadCells{i});end;