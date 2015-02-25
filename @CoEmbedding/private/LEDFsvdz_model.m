function [score,X,Y,BB, Rqz, Rq] = LEDFsvdz_model( x,R,k,obj_option, obj_para )

% x(1,2)    --- eta_1, eta_2
% x(3)      --- alpha
% x(4)      --- bta
% x(5)      --- r
% BB        --- singular values


Rq = ACAS_R(R, obj_option, obj_para);

eta1=x(1);
eta2=x(2);
zeta=x(3:k+2);
bta=x(k+3:end)';


Dr=sum(R,2);             % row sum diagonal matrix
Dc=sum(R,1);             % column sum diagonal matrix


Rx=bsxfun(@times,R,Dr.^eta1);
Pc=sum(Rx,1);

Ry=bsxfun(@times,R,Dc.^(eta2));
Qr=sum(Ry,2);

%% LED
T=bsxfun(@times,R,sqrt(Dc.^(eta2)./Pc));
[U,S,V]=svds(bsxfun(@times,T,sqrt(Dr.^(eta1)./Qr)),k+1,'L');                   % S include the nonnegative singular value in decresing order
B=diag(S).^2;


p=2:k+1 ;                  %2:k+1;
BB=B(p);                      % the k second largest eigenvalue vector



%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Xp=bsxfun(@times,sqrt(Dr.^(-eta1)./Qr),U(:,p));
gamma = sqrt(1./(diag(U(:,p)'*bsxfun(@times, Dr.^(-eta1), U(:,p)))));
Xp = bsxfun(@times, Xp, gamma');


Yp = bsxfun(@times,sqrt(Dc.^(-eta2)./Pc)',V(:,p));
    
alpha = gamma'.*zeta;
Yp = bsxfun(@times, Yp, alpha);

X=bsxfun(@times,Xp,BB.^(bta)');
Y=bsxfun(@times,Yp,BB.^(bta)');


[score, Rqz, Rq] =KNN_obj(Rq, X, Y, obj_option, obj_para);




end

