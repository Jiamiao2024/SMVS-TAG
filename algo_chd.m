function [UU,V,A,W,Z,alpha,iter,obj,res] = algo_chd(X,Y,lambda,numanchor,d,p,beta)
% m      : the number of anchor. the size of Z is m*n.
% lambda : the hyper-parameter of regularization term.
% X      : n*di

%% initialize
% maxIter = 50 ; % the number of iterations
c = length(unique(Y));
m = numanchor;
numview = length(X);
numsample = size(X{1},1);

W = cell(numview,1);            % dt * dt
A = cell(numview,1);            % dt * m
Z = zeros(m,numsample);         % m  * N
J = cell(numview,1);
res = zeros(101,1);
obj = zeros(1,51);

for t = 1:numview
    dt = size(X{t},2);
    W{t} = zeros(dt,d);
    A{t} = zeros(d,m);         % d  * m
    Y2{t} = zeros(d,m);
    J{t} = zeros(d,m);

    [X{t}, ~] = mapstd(X{t}',0,1); % turn into d*n
end
Z(:,1:m) = eye(m);
sX = [m, d, numview];
rho = 1e-4; max_rho = 10e12; pho_rho = 1.1;
tol = 1e-2;
weight_vector = ones(1,numview)';


alpha = ones(1,numview);

opt.disp = 0;
A_tensor  = zeros(m, d, numview);
J_tensor  = zeros(m, d, numview);
flag = 1;
iter = 0;

%%
while flag
    iter = iter + 1;
    A_Ten_pre =  A_tensor;
    J_Ten_pre = J_tensor;

    %% optimize W_i

    parfor iv=1:numview
        AZ{iv} = A{iv}*Z;
        C{iv} = X{iv}*AZ{iv}';
        [U,~,V] = svd(C{iv},'econ');
        W{iv} = U*V';
    end

    %% optimize A
    for ia = 1:numview
        al2 = alpha(ia)^2;
        D{ia} = al2 * beta * 2 *W{ia}' * X{ia} * Z';
        DD{ia} = -Y2{ia} - rho*J{ia};
        Dtem{ia} = D{ia}+DD{ia};
        [U,~,V] = svd(Dtem{ia},'econ');
        A{ia} = U*V';
    end

    %% optimize Z
    sumAlpha = sum(alpha.^2);
    % % QP
    % Sbar=[];
    % H = 2*sumAlpha*A'*A+2*lambda*eye(m);
    H = 2*sumAlpha*beta*eye(m)+2*lambda*eye(m);
    H = (H+H')/2;
    % [r,q] = chol(H);

    options = optimset( 'Algorithm','interior-point-convex','Display','off'); % Algorithm 默认为 interior-point-convex
    parfor ji=1:numsample
        ff=0;
        for j=1:numview

            C = W{j} * A{j};
            ff = ff - 2*beta*X{j}(:,ji)'*C;
        end
        Z(:,ji) = quadprog(H,ff',[],[],ones(1,m),1,zeros(m,1),ones(m,1),[],options);
    end


    %%  optimize J{v}
    for v =1:numview
        QQ{v}=(A{v}' - Y2{v}'/rho);  %m*d
    end
    Q_tensor = cat(3,QQ{:,:}); %m*d*V
    Qg = Q_tensor(:);
    [myj, ~] = wshrinkObj_weight_lp(Qg, weight_vector./rho,sX, 0,3,p);
    J_tensor = reshape(myj, sX);%m*d*V
    for ii=1:numview
        J{ii} = J_tensor(:,:,ii); %m*d
        J{ii} = J{ii}';
    end


    %% optimize Y and  penalty parameters
    for ii=1:numview
        Y2{ii} = Y2{ii} + rho*(J{ii}-A{ii});
    end
    rho = min(rho*pho_rho, max_rho);

    term_rec = 0;
    term_Z = 0;
    for v = 1:numview
        term_rec = term_rec + beta * norm(X{v} - W{v}*A{v}*Z, 'fro')^2;
    end
    term_Z = lambda * norm(Z, 'fro')^2;
    A_tensor = cat(3, A{:});              % d × m × V
    A_tensor = permute(A_tensor, [2,1,3]);% m × d × V
    term_A = Sp_A(A_tensor, p);

    obj(iter) = term_rec + term_Z + term_A;
    chgA = 0;% max(abs(A_Ten_pre(:)-A_tensor(:)));
    chgJ = 0;%max(abs(J_Ten_pre(:)-J_tensor(:)));
    chgA_J = max(abs(A_tensor(:)-J_tensor(:)));
    chg = max([chgA chgJ chgA_J]);

    if chg<tol || iter>50
        [UU,~,V]=svd(Z','econ');
        flag = 0;
    end

end



