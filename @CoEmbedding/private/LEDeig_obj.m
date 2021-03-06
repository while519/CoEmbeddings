function [score] = LEDeig_obj(x,R,k,obj_option, obj_para )

% x(:,1,2)    --- eta_1, eta_2
% x(:,3)      --- alpha
% x(:,4)      --- bta



Rq = ACAS_R(R, obj_option, obj_para);

eta1=x(:,1);
eta2=x(:,2);
t1=x(:,3);
bta=x(:,4:end);

len=size(x,1);

Dr=sum(R,2);             % row sum diagonal matrix              column
Dc=sum(R,1);             % column sum diagonal matrix           row

score=zeros(1,len);

parfor ii=1:len
    
    Rx=bsxfun(@times,R,Dr.^(eta1(ii)));
    Pc=sum(Rx,1);                                               %row
    
    Ry=bsxfun(@times,R,Dc.^(eta2(ii)));
    Qr=sum(Ry,2);                                               %column
    
    %% LED
    T1 = bsxfun(@rdivide, Ry, Qr);
    T2 = bsxfun(@rdivide, Rx, Pc);
    
    [U,D]=eigs(T1*T2',k+1);                   % S include the nonnegative eigen values in decresing order
    [ U, B] = EigSort( U, D, 'descend');
    
    
    p= 2:k+1;
    
    BB=B(p)';                      % the k second largest eigenvalue vector
    sigma = BB./BB(1);
    
    
    
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Xp = bsxfun(@times,sigma.^bta(ii,:),U(:,p));
    gamma = sqrt(1./(diag(U(:,p)'*bsxfun(@times, Qr, U(:,p)))));
    X = bsxfun(@times, Xp, gamma');
    
    Y = bsxfun(@rdivide,Rx'*X*t1(ii),Pc');
    Y = bsxfun(@rdivide, Y, sqrt(BB));
    
    

    
    score(ii) = KNN_obj(Rq, X, Y, obj_option, obj_para);
end

end

