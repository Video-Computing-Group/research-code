%% Matlab Toolbox for Submodular Function Optimization (v 2.0)
% Tutorial script and all implementations Andreas Krause (krausea@gmail.com)
%
% This is the Octave version.  More information can be found under
% sfo_tutorial.


do_plot = 0; %switch visualization on or off

disp(' ');
disp('---------------------------------------------------------');
disp('---------------------------------------------------------');
disp('Welcome to the submodular function optimization tutorial');
disp('---------------------------------------------------------');
disp('---------------------------------------------------------');
disp(' ')
if ~do_plot
    disp('Visualization disabled');
end
pause

%% INITIALIZATION

% We first define several examples of submodular functions
% In order to obtain more information about submodular function objects,
% type 'help sfo_fn'
%
% First some functions for experimental design on a Gaussian Process
% trained on pH data from a lake in Merced, California

% load the data: Contains covariance matrix merced_data.sigma, 
% and locations (coordinates) merced_data.coords
load merced_data;
% the ground set for the submodular functions
V_sigma = 1:size(merced_data.sigma,1); 

% Mutual information: F_mi(A) = H(V\A) - H(V\A | A)
F_mi = sfo_octavize(sfo_fn_mi(merced_data.sigma,V_sigma));
% Mutual information: F_mi(A) = H(V\A) - H(V\A | A)
F_ig = sfo_octavize(sfo_fn_infogain(merced_data.sigma,V_sigma,.01));
% Variance reduction: F_var(A) = Var(V)-Var(V | A)
F_var = sfo_octavize(sfo_fn_varred(merced_data.sigma,V_sigma));

% Helper function for evaluating the maximum remaining variance,
% eval_maxvar(A) = max_s Var(s | A);
eval_maxvar = @(A) sfo_eval_maxvar(merced_data.sigma,A); 

% *Note:* mutual information is submodular, but not monotonic! But for when
% selecting sets of size k << n (where n is the total # of elements) it's
% approximately monotonic, and that's enough (see Krause et al., JMLR '08)
% Also, variance reduction is not always submodular (see Das & Kempe, STOC '08)

% -----------------------------------------------------
% Now the 2 element example function from the tutorial slides at www.submodularity.org
% 'F{[]} = 0, F({a}) = -1, F({b}) = 2, F({a,b}) = 0';
F_ex = sfo_octavize(sfo_fn_example);
V_ex = 1:2;

% -----------------------------------------------------
% now some cut functions on directed and undirected graphs
% first define adjacency matrices
G_dir=[0 1 1.2 0 0 0; 1.3 0 1.4 0 0 0; 1.5 1.6 0 1.7 0 0; 0 0 1.8 0 1.9 2; 0 0 0 2.1 0 2.2; 0 0 0 2.3 2.4 0];
G_un=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
% here's the cut functions:
F_cut_dir = sfo_octavize(sfo_fn_cutfun(G_dir));
F_cut_un = sfo_octavize(sfo_fn_cutfun(G_un));
% and the ground set V_G
V_G = 1:6;

% -----------------------------------------------------
% now an objective function for modeling sensor detections
n_sensors = 100;
n_targets = 50;
detprob = rand(n_sensors,n_targets);
V_det = 1:n_sensors;
F_det = sfo_octavize(sfo_fn_detect(detprob,V_det));

% -----------------------------------------------------
% a test function by Iwata for testing submodular minimization algorithms
% described, e.g., in Fujishige '06
V_iw = 1:100;
F_iw = sfo_octavize(sfo_fn_iwata(length(V_iw)));

%% MINIMIZATION OF SUBMODULAR FUNCTIONS
disp(' ')
disp('---------------------------------------------------------');
disp('---------------------------------------------------------');
disp(' MINIMIZATION OF SUBMODULAR FUNCTIONS ');
disp('---------------------------------------------------------');
disp('---------------------------------------------------------');
pause

%% Queyranne's algorithm
disp(' ');
disp('---------------------------------------------------------');
disp('We will first explore Queyranne''s algorithm for minimizing a symmetric submodular function');
disp(' ');
disp('Example: F is a cut function on an undirected graph with adjacency matrix ');
G_un
disp(' ')
F = F_cut_un; V = V_G;
disp('A = sfo_queyranne(F,V)');
A = sfo_queyranne(F,V)
disp(' ')
pause

%% Minimizing general submodular functions

disp(' ');
disp('---------------------------------------------------------');
disp('Now let''s explore general submodular function minimization using the min-norm-point algorithm');
disp(' ');
disp('Example: Iwata''s test function');
F = F_iw; V = V_iw;
disp('A = sfo_min_norm_point(F,V)')
pause
A = sfo_min_norm_point(F,V)
disp(' ')
pause



%% Finding a minimum s-t-cut using general submodular function minimization

disp(' ');
disp('---------------------------------------------------------');
disp('Here''s how we can find submodular s-t-mincuts,');
disp('i.e., solutions to min_A F(A) s.t. s in A, t not in A');
disp('Example: s-t-cuts on a directed graph with adjacency matrix');
G_dir
disp(' ');
F = F_cut_dir; V = V_G;
disp('Separating node 1 from 6');
disp('A16 = sfo_s_t_mincut(F,V,1,6)');
pause
A16 = sfo_s_t_mincut(F,V,1,6)
disp('Separating node 4 from 6');
disp(' ');
disp('A46 = sfo_s_t_mincut(F,V,4,6)');
A46 = sfo_s_t_mincut(F,V,4,6)
disp(' ');
pause

%% Clustering using mutual information

randn('state',0); rand('state',0);

n = 8;

disp(' ')
disp('---------------------------------------------------------');
disp('Now we use the greedy splitting algorithm for clustering')
disp(sprintf('First generate %d data points',3*n));
disp(' ');

C1 = randn(2,n)+5;
C2 = randn(2,n);
C3 = randn(2,n)-5;

X = [C1';C2';C3'];

if do_plot
    figure
    plot(X(:,1),X(:,2),'b.'); hold on
    pause
end
    
disp('Now use a Gaussian kernel function');
disp('D = sfo_dist(X)/2; sigma_cl = exp(-D.^2)+eye(3*n)*.01;');
D = sfo_dist(X)/2;
sigma_cl = exp(-D.^2)+eye(3*n)*.01;

disp('Use entropy as energy function');
disp('E = sfo_octavize(sfo_fn_entropy(sigma_cl,A));');
V = 1:(3*n); 
E = sfo_octavize(sfo_fn_entropy(sigma_cl,A));
disp(' ');
disp('Now do greedy splitting with 2 clusters')
disp('P = sfo_greedy_splitting(E,V,2)');
pause
P = sfo_greedy_splitting(E,V,2)

disp('Now plot clusters');

if do_plot
    clf
    plot(X(P{1},1),X(P{1},2),'bx'); hold on
    plot(X(P{2},1),X(P{2},2),'ro');
    disp(' ')
    pause
end

disp('Now do greedy splitting with 3 clusters')
disp('P = sfo_greedy_splitting(E,V,3)');
pause
P = sfo_greedy_splitting(E,V,3)

disp('Now plot clusters');

if do_plot
    clf
    plot(X(P{1},1),X(P{1},2),'bx'); hold on
    plot(X(P{2},1),X(P{2},2),'ro');
    plot(X(P{3},1),X(P{3},2),'ks'); hold off
    disp(' ')
end
pause

%% Image denoising
randn('state',0); rand('state',0);

disp(' ');
disp('---------------------------------------------------------');
disp('We will now do inference in a Markov Random Field using ');
disp('submodular function minimization');
disp(' ');
D = 40; S = 10; noiseProb = 0.2;
img = zeros(D);
img(S:(D-S),S:(D-S))=1;
mask = rand(size(img))>(1-noiseProb);
imgN = img; imgN(mask)=1-imgN(mask);
V = 1:numel(img);

if do_plot
    figure; 
    subplot(131); imshow(img); title('original image')
    subplot(132); imshow(imgN); title('noisy image')
end

% potentials for use in the ising model
coeffPix = 1; coeffH = 1; coeffV = 1; coeffD = .0;

disp('now we define a submodular function on the noisy image imgN');
F = sfo_octavize(sfo_fn_ising(imgN, coeffPix, coeffH, coeffV, coeffD));
disp('F = sfo_octavize(sfo_fn_ising(imgN,A, coeffPix, coeffH, coeffV, coeffD))');

disp(' ');
disp('now minimize using min norm point algorithm')
if do_plot
    pause
    subplot(133)
    callback = @(sol) sfo_min_norm_point_tutorial_callback_image(sol,zeros(size(imgN)));
else
    callback = [];
end
tic 
Ainit = find(imgN(:)); %use for initialization
opt = sfo_opt({'minnorm_init',Ainit,'minnorm_stopping_thresh',3.999,'minnorm_callback',callback});
AD = sfo_min_norm_point(F,V,opt); %allow to be at most 1 pixels off
toc
if do_plot
    imgD = zeros(size(img)); imgD(AD)=1;
    imshow(imgD)
    title('reconstructed image')
end
disp(' ')
pause




%% MAXIMIZATION OF SUBMODULAR FUNCTIONS
disp(' ')
disp('---------------------------------------------------------');
disp('---------------------------------------------------------');
disp(' MAXIMIZATION OF SUBMODULAR FUNCTIONS ');
disp('---------------------------------------------------------');
disp('---------------------------------------------------------');
pause

%% The lazy greedy algorithm

disp(' ');
disp('---------------------------------------------------------');
disp('We will now explore the lazy greedy algorithm');
disp(' ');
k = 15;
disp('Example: Experimental Design. Want to predict pH values on Lake Merced');
disp(sprintf('We greedily pick %d sensor locations for maximizing mutual information',k));
F = F_mi; V = V_sigma;
disp('[A,scores,evals] = sfo_greedy_lazy(F,V,k);');
pause
[A,scores,evals] = sfo_greedy_lazy(F,V,k);
A
nevals = length(V):-1:(length(V)-k+1);
disp(sprintf('Lazy evaluations: %d, naive evaluations: %d, savings: %f%%',sum(evals),sum(nevals),100*(1-sum(evals)/sum(nevals))));
disp(' ')

if do_plot
    figure
    plot(merced_data.coords(:,1),merced_data.coords(:,2),'k.'); hold on
    plot(merced_data.coords(A,1),merced_data.coords(A,2),'bs','markerfacecolor','blue');
    xlabel('Horizontal location along transect'); ylabel('Vertical location (depth)'); title('Chosen sensing locations (blue squares)');
    pause
end

disp('Let''s now compute online bounds.  These are data dependent bounds');
disp('that can be computed after running the greedy algorithm (and are often');
disp('tighter than the (1-1/e) "offline" bound).');
bound = sfo_maxbound(F,V,A,k);
disp(sprintf('Greedy score F(A) = %f; \nNemhauser (1-1/e) bound: %f; \nOnline bound: %f',scores(end),scores(end)/(1-1/exp(1)),bound));
pause

if do_plot
    disp('will now plot the greedy scores');
    figure
    plot(0:k,[0 scores]); xlabel('Number of elements'); ylabel('submodular utility');
    pause
end

disp(' ');
disp('here''s how we can do greedy optimization over a matroid')
disp('define a function that takes a set A and outputs 1 if A is independent, 0 otherwise')
opt = sfo_opt({'greedy_check_indep',(@(A) (length(A)<=k))});
disp('Example: Uniform matroid: A independent <=> length(A)<=k')
disp('opt = sfo_opt({''greedy_check_indep'',@(A) (length(A)<=k)})');
disp(' ');
disp('Now let''s run the greedy algorithm, given infinite budget:');
disp('[A,scores,evals] = sfo_greedy_lazy(F,V,inf,opt)');
pause
[A,scores,evals] = sfo_greedy_lazy(F,V,inf,opt);
A


%% The lazy greedy coverage algorithm

disp(' ')
disp('---------------------------------------------------------');
disp('Now for submodular coverage');
disp('Example: Picking best sensor locations to achieve a specified amount of mutual information');
disp(' ');
F = F_mi; V = V_sigma;
C = 1:(length(V));%make up some cost function
opt = sfo_opt({'cost',C}); 
Q = 5;
disp(sprintf('Covering quota Q = %f',Q));
disp('[A,stat] = sfo_cover(F,V,Q,opt)');
pause
[A,stat] = sfo_cover(F,V,Q,opt);
disp(sprintf('Coverage possible: %d, Cost: %f',stat,sum(C(A))));
A
pause
Q = 30;
disp(sprintf('Covering quota Q = %f',Q));
disp('[A,stat] = sfo_cover(F,V,Q,opt)');
[A,stat] = sfo_cover(F,V,Q,opt);
disp(sprintf('Coverage possible: %d, Cost: %f',stat,length(A)));
A
disp(' ');
pause

%% The CELF algorithm for budgeted maximization
disp(' ');
disp('---------------------------------------------------------');
disp('Now we test the CELF algorithm');
disp(' ');
F = F_mi; V = V_sigma; 
% make up some cost function
C = 1:length(V_sigma);
opt = sfo_opt({'cost',C}); 

disp('First use a small budget, B = 2. Here, unit cost is better');
sfo_celf(F,V,2,opt)
pause
disp('First use a small budget, B = 15. Here, cost/benefit is better');
sfo_celf(F,V,15,opt)
disp(' ')
pause

%% The local search algorithm for nonnegative maximization
disp(' ');
disp('---------------------------------------------------------');
disp('Now we test the Local Search algorithm');
disp(' ');
F_cost = sfo_fn_wrapper(@(A) length(sfo_unique_fast(A)));
F = sfo_octavize(sfo_fn_lincomb({F_var,F_cost},[1 -.001]));
V = V_sigma; 

disp('Optimizing variance reduction minus sensing cost')
fprintf('For this problem instance, F(V)=%f > 0\n',F(V))
disp(' ')
disp('Let''s run local search: ')
disp('sfo_ls_lazy(F,V)')
sfo_ls_lazy(F,V)
disp(' ')
pause

%% Welfare and balance problems
disp(' ')
disp('---------------------------------------------------------');
disp('We will now compute balanced sensor placements');
V = V_sigma;
F = F_var;
disp('pick 20 sensor locations to activate at 5 timeslices')
disp('optimize the sum of the performance:')
Fs = {sfo_octavize(sfo_fn_varred(merced_data.sigma,V_sigma)),sfo_octavize(sfo_fn_varred(merced_data.sigma,V_sigma)),sfo_octavize(sfo_fn_varred(merced_data.sigma,V_sigma)),sfo_octavize(sfo_fn_varred(merced_data.sigma,V_sigma)),sfo_octavize(sfo_fn_varred(merced_data.sigma,V_sigma))};
[A_sum,scores_sum] = sfo_greedy_welfare(Fs,V,20);
pause
disp('now optimize the min of the performance:')
[A_min,scores_min] = sfo_balance(F,V,5,20);
disp(sprintf('scores sum = %s',sprintf('%f ',scores_sum)));
disp(sprintf('scores min = %s',sprintf('%f ',scores_min)));

pause

%% Submodular-supermodular procedure of Narasimhan & Bilmes
randn('state',0); rand('state',0);
disp(' ');
disp('---------------------------------------------------------');
disp('We will now use the submodular-supermodular procedure to do experimental design');
disp(' ');
disp('Want to choose subset A of sensor locations in the Lake Merced ph-estimation data set')
disp('that maximizes MI(A)-|A|, where MI is the mutual information criterion');
disp('We set F(A) = |A|, and G(A) = MI(A), and minimize F(A)-G(A)');
disp('G = F_mi; F = sfo_fn_wrapper(@(A) length(A))');
V = V_sigma;
G = F_mi;
F = sfo_fn_wrapper(@(A) length(sfo_unique_fast(A)));
disp(' ');
disp('A = sfo_ssp(F,G,V)')
pause
A = sfo_ssp(F,G,V_sigma)
disp(sprintf('\n\nMutual information MI(A) = %f, Cost |A| = %d',G(A),F(A)));
disp(' ');

if do_plot
    figure
    plot(merced_data.coords(:,1),merced_data.coords(:,2),'k.'); hold on
    plot(merced_data.coords(A,1),merced_data.coords(A,2),'b*','markerfacecolor','blue','markersize',10);
    xlabel('Horizontal location along transect'); ylabel('Vertical location (depth)'); title('Chosen locations by SSSP (blue stars)');
end
pause

%% The Data-Correcting algorithm for maximizing general submodular functions

disp(' ');
disp('---------------------------------------------------------');
disp('We will now use the Data Correcting algorithm for finding');
disp('a maximum directed cut in a graph');
disp(' ');
F = F_cut_dir; V = V_G;
disp('A = sfo_max_dca_lazy(F,V)');
pause
A = sfo_max_dca_lazy(F,V);
disp(sprintf('Cut value = %f',F(A)));
disp(' ');
pause
disp('Now let''s find the best cut of size 1:');
disp('Define a new submodular function that''s -inf if |A|>k:');
disp('FT = sfo_fn_wrapper(@(A) F_cut_dir(A) - 1e10*max(0,length(A)-1))')
FT = sfo_fn_wrapper(@(A) F_cut_dir(A) - 1e10*max(0,length(A)-1));
disp('A = sfo_max_dca_lazy(FT,V)');
pause
A = sfo_max_dca_lazy(FT,V);
disp(sprintf('Cut value = %f',FT(A)));
disp(' ');

pause

%% The Saturate algorithm for robust optimization

disp(' ')
disp('---------------------------------------------------------');
disp('Now let''s use the Saturate algorithm to minimize ')
disp('worst-case variance in GP regression on Merced Lake')
disp(' ');
V = V_sigma;
k = 10;
disp('Run greedy algorithm: ')
disp('AG = sfo_greedy_lazy(F_var,V,k)');
AG = sfo_greedy_lazy(F_var,V,k)
disp(' ');
if do_plot
    figure
    plot(merced_data.coords(:,1),merced_data.coords(:,2),'k.'); hold on
    plot(merced_data.coords(AG,1),merced_data.coords(AG,2),'bs','markerfacecolor','blue');
    xlabel('Horizontal location along transect'); ylabel('Vertical location (depth)'); title('Greedy locations (blue squares) and Saturate locations (green diamonds)');

    pause
end

disp('Run Saturate algorithm: ')
disp('AS = sfo_saturate(F_var,V,k,''min_thresh'')');
pause
AS = sfo_saturate(sfo_fn_varred(merced_data.sigma,V_sigma),V,k,'minthresh');
if do_plot
    plot(merced_data.coords(AS,1),merced_data.coords(AS,2),'gd','markerfacecolor','green');hold off
end
disp(' ');
disp(sprintf('max remaining variance: Greedy = %f, Saturate = %f',eval_maxvar(AG),eval_maxvar(AS)));
pause


%% The pSPIEL Algorithm for trading off accuracy and communication cost

disp(' ')
disp('---------------------------------------------------------');
disp('Now let''s use the pSPIEL algorithm to trade off informativeness and')
disp('communication cost in GP regression on Merced Lake')
disp(' ');
V = V_sigma;
D = merced_data.dists + ones(length(V)); %cost of links + cost of nodes
Q = 0.6*F_var(V); %want 80% of optimal variance reduction

disp('Run greedy algorithm: ')
disp('AG = sfo_cover(F_var,V,Q)');
AG = sfo_cover(F_var,V,Q)

disp('AP = sfo_pspiel(F_var,V,Q,D)');
AP = sfo_pspiel(F_var,V,Q,D)

utility_greedy = F_var(AG);
utility_pspiel = F_var(AP);
[cost_greedy edges_greedy steiner_greedy] = sfo_pspiel_get_cost(AG,D);
[cost_pspiel edges_pspiel steiner_pspiel]= sfo_pspiel_get_cost(AP,D);

if do_plot
    subplot(211)
    sfo_plot_subgraph(merced_data.coords,edges_greedy,steiner_greedy)
    title('Greedy-connect');
    subplot(212)
    sfo_plot_subgraph(merced_data.coords,edges_pspiel,steiner_pspiel)
    title('pSPIEL');
end

disp(sprintf('Greedy-connect: Utility = %f, Cost = %f. pSPIEL: Utility = %f, Cost = %f.',utility_greedy,cost_greedy,utility_pspiel,cost_pspiel));
pause

%% MISCELLANEOUS OTHER FUNCTIONS
disp(' ')
disp('---------------------------------------------------------');
disp('---------------------------------------------------------');
disp(' MISCELLANEA ');
disp('---------------------------------------------------------');
disp('---------------------------------------------------------');
pause

%% The Lovasz extension
F = F_ex; V = V_ex;
A = 2;

disp(' ')
disp('---------------------------------------------------------');
disp('Here''s how we can compute the Lovasz extension for the ');
disp('submodular function from the tutorial,');
disp('F{[]} = 0, F({a}) = -1, F({b}) = 2, F({a,b}) = 0');
disp(' ');
disp('Let''s compute the characteristic vector for set A={b}')
disp('w = sfo_charvector(V,A) ');
w = sfo_charvector(V,A)

disp(' ');
disp('Defining the Lovasz extension ');
disp('g = @(w) sfo_lovaszext(F,V,w);');
g = @(w) sfo_lovaszext(F,V,w);
disp(' ')
disp('comparing the function F(A) with Lovasz extension g(w). Should be equal');
disp(sprintf('A = {b};   F(A) = %f;   g(wA) = %f',F(A),g(w)));
disp(' ')
pause

disp('now plot the Lovasz extension');
if do_plot
    [X,Y] = meshgrid(0:.05:1,0:.05:1);
    Z = zeros(size(X));
    for i = 1:size(X,1)
        for j = 1:size(X,2)
           Z(i,j) = g([X(i,j),Y(i,j)]);
        end
    end
    figure
    surf(X,Y,Z)
    xlabel('w_{a}')
    ylabel('w_{b}')
    zlabel('g(w)')
    title('Lovasz extension for example from tutorial: F(\{\})=0,F(\{a\})=-1,F(\{b\})=2,F(\{a,b\})=0 ');
end
disp(' ')
pause

%% Bounds on optimal solution for minimization

disp(' ');
disp('---------------------------------------------------------');
disp('Here''s how we can compute bounds on the minimum solution')
disp(' ');
F = F_ex; V = V_ex;
disp('Set A={a}')
A = 1;

disp('Compute bound:')
disp('bound = sfo_minbound(F,V,A);')
bound = sfo_minbound(F,V,A);


% computing a lower bound on min F(A)
disp(sprintf('bound = %f <= min F(A) <= F(B) = %f; ==> A is optimal!',bound,F(A)));
disp(' ');
pause

%% The polyhedron greedy algorithm 
disp(' ');
disp('---------------------------------------------------------');
disp('Here''s the polyhedron greedy algorithm')
disp(' ');
F = F_ex; V = V_ex;
disp('Set A = {a}');
A = 1;

w = sfo_charvector(V,A); 
disp('xw = sfo_polyhedrongreedy(F,V,w)');
xw = sfo_polyhedrongreedy(F,V,w)
disp(' ');
pause
