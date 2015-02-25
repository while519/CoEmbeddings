
function R = CalRelationXY(X, Y, varargin);

% Computation of (dis)simialtiy features
% T. Mu, Feb. 2011

if length(varargin) ==0
    %default values
    type   = 'cosine';
    para   = [];
  elseif length(varargin) == 1
    type = varargin{1};
    para = [];
  elseif length(varargin) == 2
    type = varargin{1};
    para = varargin{2};
  end    
  

  switch lower(type)    
    case 'dot'
      R = X*Y';
    case 'cosine'
      R =X*Y';
      ind = find(R==0);
      D1 = 1./sqrt(diag(X*X'));
      D2 = 1./sqrt(diag(Y*Y'));
      R = bsxfun(@times, R, D1);
      R = bsxfun(@times, R, D2');
      R(ind) = 0;
    case {'euclidean','sqeuclidean'} 
      DOT=X*Y';
      T1=sparse(diag(diag(X*X')));
      T2=sparse(diag(diag(Y*Y')));
      R = real(sqrt(T1*ones(size(DOT))+ones(size(DOT))*T2-2*DOT));  
    case {'sqeuclidean'} 
      DOT=X*Y';
      T1=sparse(diag(diag(X*X')));
      T2=sparse(diag(diag(Y*Y')));
      R = T1*ones(size(DOT))+ones(size(DOT))*T2-2*DOT;
    case 'euclidean_sim' 
      DOT=X*Y';
      T1=sparse(diag(diag(X*X')));
      T2=sparse(diag(diag(Y*Y')));
      D = real(sqrt(T1*ones(size(DOT))+ones(size(DOT))*T2-2*DOT));
      R= max(max(D))-D;
    case 'tanimoto' % Tanimoto similarity coefficient
      R  = X * Y';
      D1 = diag(X*X');
      D2 = diag(Y*Y');
      R1 = bsxfun(@plus, -R, D1);
      R1 = bsxfun(@plus, R1, D2');
      R  = R./R1;
    case 'gaussian'
      DOT=X*Y';
      T1=sparse(diag(diag(X*X')));
      T2=sparse(diag(diag(Y*Y')));
      D = T1*ones(size(DOT))+ones(size(DOT))*T2-2*DOT;
      if isempty(para),                p = median( D(:) ) * 0.5; R = exp(-D/ (2*p^2)); 
      elseif isequal(para, 'max'),     R = exp(-D/ max(D(:)) ); 
      elseif isequal(para, 'mean'),    R = exp(-D/ mean(D(:)) );
      else,                            R = exp(-D/(2*para^2));
      end
    case 'gaussian2'
      DOT=X*Y';
      T1=sparse(diag(diag(X*X')));
      T2=sparse(diag(diag(Y*Y')));
      D = T1*ones(size(DOT))+ones(size(DOT))*T2-2*DOT;
      R  = exp(-para*D);
    case 'poly'
      DOT = X*Y';
      if isempty(para),              R = (1+DOT).^2;
      elseif length(para)==1,        R = (1+DOT).^para;  
      elseif length(para)==2,        R = (para(1)+DOT).^para(2);  
      end
    case 'lsg' %local scaling gaussian
      DOT=X*Y';
      T1=sparse(diag(diag(X*X')));
      T2=sparse(diag(diag(Y*Y')));
      D = T1*ones(size(DOT))+ones(size(DOT))*T2-2*DOT;    % Square of Euclidean distance
      [W1,W2]    = NeighborFinder(D, para); 
      Dknn1  = 1./sqrt(max(D.* W1)); Dknn2  = 1./sqrt(max(D.* W2,[],2)); 
      R = bsxfun(@times, D, Dknn1);
      R = bsxfun(@times, R, Dknn2);
      R = exp(-R); 

    case 'alw' % alternative weights
      DOT=X*Y';
      D1 = diag(X*X');
      D2 = diag(Y*Y');
      T1 = sparse(diag(D1));
      T2 = sparse(diag(D2));
      D  = T1*ones(size(DOT))+ones(size(DOT))*T2-2*DOT; 
      R  = bsxfun(@plus, zeros(size(D)), D1);
      R  = bsxfun(@plus, R, D2'); 
      R = 1./(D./R + para);
      

    case 'correlation'
      d = size(X,2);
      R =X*(sparse(eye(d)) - 1/d)*Y';
      D1 = 1./sqrt(diag(X*(sparse(eye(d)) - 1/d)*X'));
      D2 = 1./sqrt(diag(Y*(sparse(eye(d)) - 1/d)*Y'));
      R = bsxfun(@times, R, D1);
      R = bsxfun(@times, R, D2');  
  end
  
end