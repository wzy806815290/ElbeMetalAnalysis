%% Example of the source apportionmnet in the Elbe River downstream
%% Initialize the factors of q, J, T, and muP
q = 4;
readPath = strcat("PMF_D.xls");
[num, txt, all] = xlsread(readPath, 'conc');
info = txt(2:end, 1:2);
num_available = num(~any(num==-999,2),:);
info_available = info(~any(num==-999,2),:);

J = size(num_available, 2); % number of variables
T = size(num_available, 1); % number of samples

muPRaw = [0	0.000172318	2.33207E-05	5.32465E-05	0.863645591	0.000864468	0.127921019	0.000886799	1.78366E-07	0.006433059;
0.008198366	0.032747196	0.000534554	0	0	0.008852069	0.83603011	0	0.00083863	0.112799076;
0.051779395	0	0.002304526	0.018712084	0	0.089163786	0.342274348	0.073835949	0.000146848	0.421783064;
0.001749356	0.005124808	0	0.002453168	0.950963792	1.10463E-05	0.033702857	2.52677E-05	1.11431E-07	0.005969594];
%muPRaw = ones(q, J) * 0.5;
a = 0;
b = 0;
c = 0;
d = 0;
muP = muPRaw;

%% Run the BMRM model without pre-assignment of zero elements in F (the base model)
% BNFA(Y, q, muP, nBurnIn, nIter)
% Y: T by J data matrix (T: number of observations, J: number of variables)
% q: number of factors (number of major pollution sources)
% muP: 	prior mean matrix (of size q by J) for the source composition matrix P. 
%       Zeros need to be assigned to prespecified elements of muP to satisfy the identifiability condition C1. 
%       For the remaining free elements, any nonnegative numbers (between 0 and 1 preferably) can be assigned. 
%       If no or an insufficient number of zeros are preassigned in muP, estimation can still be performed but the resulting estimates may be subject to rotational ambiguity.
% nBurnIn: number of iterations for the burn-in period in MCMC
% nIter: number of iterations in the sampling phase of MCMC
[Phat,Ahat,Sigmahat,stdP,stdA,stdSigma,PSnor,ASnor,SigmaS] = BNFA(num_available,q,muP,10000,10000); 
for i = 1:size(PSnor, 3)
    Ds(i) = -2*sum(log(normpdf(0,0,SigmaS(i, :))));
end
D_mean = -2*sum(log(normpdf(0,0,Sigmahat)));
DIC = 2*mean(Ds)-D_mean;
BIC = D_mean + (J * (q + 2) - q * (q - 1)) * log(T);
QLoss = sum(sum((num_available-Ahat * Phat).^2));
savePath = strcat("D", num2str(q), "_ndist", num2str(a), num2str(b), num2str(c), num2str(d), "_result.mat");
save(savePath, 'q', 'muP', 'Phat', 'Ahat', 'Sigmahat', 'stdP', 'stdA', 'stdSigma', 'DIC', 'BIC', 'QLoss');

%% Run the BMRM model with pre-assignment of zero elements in F (the candidate model)
x1 = combntns([1, 2, 3, 4, 6, 8, 9], q - 1);
x2 = combntns([3, 4, 5, 8, 9], q - 1);
x3 = [2, 5, 9];
x4 = combntns([1, 3, 6, 8, 9], q - 1);
for a = 1:size(x1, 1)
    for b = 1:size(x2, 1)
        for c = 1:size(x3, 1)
            for d = 1:size(x4, 1)
                muP = muPRaw;
                muP(1, x1(a, :)) = 0;
                muP(2, x2(b, :)) = 0;
                muP(3, x3(c, :)) = 0;
                muP(4, x4(d, :)) = 0;

                [Phat,Ahat,Sigmahat,stdP,stdA,stdSigma,PSnor,ASnor,SigmaS] = BNFA(num_available,q,muP,10000,10000);
                for i = 1:size(PSnor, 3)
                    Ds(i) = -2*sum(log(normpdf(0,0,SigmaS(i, :))));
                end
                D_mean = -2*sum(log(normpdf(0,0,Sigmahat)));
                DIC = 2*mean(Ds)-D_mean;
                BIC = D_mean + (J * (q + 2) - q * (q - 1)) * log(T);
                QLoss = sum(sum((num_available-Ahat * Phat).^2));
                savePath = strcat("D", num2str(q), "_ndist", num2str(a), num2str(b), num2str(c), num2str(d), "_result.mat");
                save(savePath, 'q', 'muP', 'Phat', 'Ahat', 'Sigmahat', 'stdP', 'stdA', 'stdSigma', 'DIC', 'BIC', 'QLoss');
            end
        end
    end
end



