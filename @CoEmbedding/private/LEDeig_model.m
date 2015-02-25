function [score,X,Y,BB, Rqz, Rq] = LEDeig_model( x,R,k,obj_option, obj_para )

% x(1,2)    --- eta_1, eta_2
% x(3)      --- alpha
% x(4)      --- bta
% x(5)      --- r
% BB        --- eigen values


Rq = ACAS_R(R, obj_option, obj_para);

eta1=x(1);
eta2=x(2);
t1=x(3);
bta=x(4:end)';


Dr=sum(R,2);             % row sum diagonal matrix     column vector
Dc=sum(R,1);             % column sum diagonal matrix  row vector


Rx=bsxfun(@times,R,Dr.^eta1);
Pc=sum(Rx,1);

Ry=bsxfun(@times,R,Dc.^(eta2));
Qr=sum(Ry,2);


%% LED
    T1 = bsxfun(@rdivide, Ry, Qr);
    T2 = bsxfun(@rdivide, Rx, Pc);
    
    
    [U,D]=eigs(T1*T2',k+1);                   % S include the nonnegative eigen values in decresing order
    [ U, B] = EigSort( U, D, 'descend');
    
    
    p= 2:k+1;
    
    BB=B(p);                      % the k second largest eigenvalue vector
    sigma = BB/BB(1);
    
    
    
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Xp = bsxfun(@times,(sigma.^bta)',U(:,p));
    gamma = sqrt(1./(diag(U(:,p)'*bsxfun(@times, Qr, U(:,p)))));
    X = bsxfun(@times, Xp, gamma');
    
    Yp = bsxfun(@rdivide,Rx'*X*t1,Pc');
    Y = bsxfun(@rdivide, Yp, sqrt(BB'));

    
[score, Rqz, Rq] = KNN_obj(Rq, X, Y, obj_option, obj_para);




end

