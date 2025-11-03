function[U] = RSTR(X, v, n, dd, c,alpha, beta,gamma,omega,scalar,dataset_name)
MaxIter = 15;
U = cell(1,v);
Z = cell(1,v);
E = cell(1,v);
F = cell(1,v);
Lv=cell(1,v);  
Q=cell(1,v);  
Dv = cell(1,v);  
for  i = 1:v
    Q{i} = eye(dd(i));
    U{i} = ones(dd(i),c);
    Z{i} = ones(n,n);
    E{i} = zeros(n,n);
    Dv{i} = eye(n);
end
options = [];
options.NeighborMode = 'KNN';
options.k = 5;
options.WeightMode = 'Binary';   
for i = 1:v
    TZ = constructW(X{i}'*U{i},options);
    TZ = full(TZ);
    Z1 = TZ-diag(diag(TZ));
    TZ = (Z1+Z1')/2;
    try
        opts.tol = 1e-4;
        [F{i},~] = eigs(TZ,c,'sa',opts);
    catch ME
        if (strcmpi(ME.identifier,'MATLAB:eig:NoConvergence'))
            opts.tol = 1e-4;
            [F{i},~] = eigs(Lv, eye(size(Lv)),num_cluster,'sa',opts);
        else
            rethrow(ME);
        end
    end
end
J = cell(1,v);
O = cell(1,v);
for i = 1:v
    J{i} = zeros(n,n);
    O{i} = zeros(n,n);  
end
sX = [n, n, v];
pho = 2;
mu = 0.1;
for iter = 1:MaxIter
    for iterv = 1:v
        LG = (eye(n) - (Z{iterv}+E{iterv}));
        LG = LG * LG';
        LG = (LG + LG') / 2;
        [Y, ~, ~]=eig1(LG, c, 0);
        U{iterv}=(X{iterv}*X{iterv}'+alpha*Q{iterv})\(X{iterv}*Y);
        Ui=sqrt(sum(U{iterv}.*U{iterv},2)+eps);
        q=0.5./Ui;
        Q{iterv}=diag(q);
    end
    for iterv = 1:v
        tmepe1 = 2*beta*Dv{iterv}+2*X{iterv}'*U{iterv}*U{iterv}'*X{iterv};
        tmepe2 = 2*X{iterv}'*U{iterv}*U{iterv}'*X{iterv}-2*X{iterv}'*U{iterv}*U{iterv}'*X{iterv}*Z{iterv};
        tempe  = tmepe1\tmepe2;
        E{iterv}=tempe;
        Ei=sqrt(sum(E{iterv}.*E{iterv},2)+eps);
        e=0.5./Ei;
        Dv{iterv}=diag(e);
    end
    for iterv = 1:v 
        M = pinv(2*X{iterv}'*U{iterv}*U{iterv}'*X{iterv}+mu*eye(n)) * (2*X{iterv}'*U{iterv}*U{iterv}'*X{iterv} -2*X{iterv}'*U{iterv}*U{iterv}'*X{iterv}*E{iterv}- gamma*L2_distance_1(F{iterv}',F{iterv}') +mu*J{iterv}-O{iterv}) / scalar;
        if(strcmp(dataset_name, 'COIL20'))
            M = double(M);
        end
        M = SimplexProj(M');
        Z{iterv} = scalar*M';
    end
    for iterv = 1:v
        sum_U = (Z{iterv}+Z{iterv}')*0.5;
        Lv = diag(sum(sum_U))-sum_U;
        Lv = (Lv+Lv')*0.5;
        try
            opts.tol = 1e-4;
            [F{iterv},~] = eigs(Lv,c,'sa',opts);
        catch ME
            if (strcmpi(ME.identifier,'MATLAB:eig:NoConvergence'))
                opts.tol = 1e-4;
                [F{iterv},~] = eigs(Lv, eye(size(Lv)),num_cluster,'sa',opts);
            else
                rethrow(ME);
            end
        end
    end
    Z_tensor = cat(3, Z{:,:});
    O_tensor = cat(3, O{:,:});
    z = Z_tensor(:);
    o = O_tensor(:);
    [j, ~] = wshrinkObj(z+1/mu*o,1/mu,sX,0,3,omega);
    J_tensor = reshape(j, sX);
    for i=1:v
        J{i} = J_tensor(:,:,i);
    end
    for i=1:v
        O{i} = O{i} + mu*(Z{i}-J{i});
    end
    mu = pho*mu;
    for iv = 1:v
        leq{iv} = Z{iv}-J{iv};
    end
    leqm = cat(3,leq{:,:});
    leqm2 = max(abs(leqm(:)));
    matchERR = max(leqm2);
end


