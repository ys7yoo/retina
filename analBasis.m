
klength = 9


subplot(211)
nkbasis = 4

ktbasprs.neye = 0; % number of "identity" basis vectors
ktbasprs.ncos = nkbasis; % Number of raised-cosine vectors to use
ktbasprs.kpeaks = [0 klength*(1 - 1.5/nkbasis)];  % Position of 1st and last bump
ktbasprs.b = 10; % Offset for nonlinear scaling (larger -> more linear)
[~,ktbasis] = makeBasis_StimKernel(ktbasprs,klength);
    
plot(ktbasis)
box off


subplot(212)
nkbasis = 8

ktbasprs.neye = 0; % number of "identity" basis vectors
ktbasprs.ncos = nkbasis; % Number of raised-cosine vectors to use
ktbasprs.kpeaks = [0 klength*(1 - 1.5/nkbasis)];  % Position of 1st and last bump
ktbasprs.b = 10; % Offset for nonlinear scaling (larger -> more linear)
[~,ktbasis] = makeBasis_StimKernel(ktbasprs,klength);
    
plot(ktbasis)
xlabel('time bin')
box off