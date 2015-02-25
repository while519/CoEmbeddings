function [score] = LEDsvdz_obj(x,R,k,obj_option, obj_para )

% x(:,1,2)    --- eta_1, eta_2
% x(:,3)      --- zeta
% x(:,4)      --- bta


Rq = ACAS_R(R, obj_option, obj_para);

eta1=x(:,1);
eta2=x(:,2);
zeta=x(:,3);
bta=x(:,4:end);

len=size(x,1);

Dr=sum(R,2);             % row sum diagonal matrix
Dc=sum(R,1);             % column sum diagonal matrix

score=zeros(1,len);

parfor ii=1:len
    
    Rx=bsxfun(@times,R,Dr.^(eta1(ii)));
    Pc=sum(Rx,1);
    
    Ry=bsxfun(@times,R,Dc.^(eta2(ii)));
    Qr=sum(Ry,2);
    
    %% LED
    T=bsxfun(@times,R,sqrt(Dc.^(eta2(ii))./Pc));
    [U,S,V]=svds(bsxfun(@times,T,sqrt(Dr.^(eta1(ii))./Qr)),k+1,'L');                   % S include the nonnegative singular value in decresing order
    B=diag(S).^2;
    
    
    p= 2:k+1;
    BB=B(p);                      % the k second largest eigenvalue vector
    
    
    
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Xp = bsxfun(@times,sqrt(Dr.^(-eta1(ii))./Qr),U(:,p));
    gamma = sqrt(1./(diag(U(:,p)'*bsxfun(@times, Dr.^(-eta1(ii)), U(:,p)))));
    Xp = bsxfun(@times, Xp, gamma');
    
    Yp = bsxfun(@times,sqrt(Dc.^(-eta2(ii))./Pc)',V(:,p));
   alpha = zeta(ii)*gamma;
    Yp = bsxfun(@times, Yp, alpha');
        
    X = bsxfun(@times, Xp, BB.^(bta(ii,:)')');
    Y = bsxfun(@times, Yp, BB.^(bta(ii,:)')');
   
    score(ii) = KNN_obj(Rq, X, Y, obj_option, obj_para);
    
end




end

