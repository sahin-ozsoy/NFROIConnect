function W = lcmvFast(data,L,fReset,para)

% 2017 Stefan Haufe

persistent Lvox;
persistent repE2;

data = data';
[M_,N,NDUM] = size(L);

out.H = eye(M_);
out.noisecov = eye(M_);
out.spatfilt = eye(M_);

if isequal(data, data')
    out.datacov = data;    
else
    if isnan(para.alpha)
      [out.datacov out.alpha] = shrinkage(data, struct('Gamma', 'auto'));
    else
      out.datacov = shrinkage(data, struct('Gamma', para.alpha));
      out.alpha = para.alpha;
    end
    data = out.spatfilt*data;
end  
    
M = size(out.spatfilt, 1);

CI = pinv(out.datacov);

if para.onedim
  W = zeros(M, N);  
else
  W = zeros(M, N, 3);  
end

if isempty(Lvox) || fReset
    out.L = L;
    for ivox = 1:N
        Lvox(:,:,ivox) = reshape(out.L(:, ivox, :), M, NDUM);
        E2 = inv(Lvox(:,:,ivox)'*Lvox(:,:,ivox));
        repE2(:,:,ivox) = repmat(sqrt(diag(E2))', M, 1);
       
    end
end

for ivox = 1:N

    % Lvox = reshape(out.L(:, ivox, :), M, NDUM);
    
    % E1 = inv(Lvox(:,:,ivox)'*CI*Lvox(:,:,ivox));
    % E2 = inv(Lvox'*CInoise*Lvox);
    
    W_ = ((Lvox(:,:,ivox)'*CI*Lvox(:,:,ivox))\Lvox(:,:,ivox)'*CI)';
    % out.W(:, ivox, :) = W_;

    W(:, ivox, :) = W_ ./ repE2(:,:,ivox);      
end

