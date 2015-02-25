function [score] = LEDFsvd_obj(x,R,k,obj_option, obj_para )

% x(:,1,2)    --- eta_1, eta_2
% x(:,3)      --- alpha
% x(:,4)      --- bta


% for ii = 1 : 2 : length(obj_para)-1
%     switch lower(obj_para{ii})
%         case 'knn'
%             knn = obj_para{ii+1};
%         otherwise
%             error( [ 'unknown input argument: ' I{ii} ] );
%     end
% end

Rq = ACAS_R(R, obj_option, obj_para);

eta1=x(:,1);
eta2=x(:,2);
t1=x(:,3:k+2);           % alpha
bta=x(:,k+3:end);

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
    Yp = bsxfun(@times,sqrt(Dc.^(-eta2(ii))./Pc)',V(:,p)*diag(t1(ii,:)));
    
    
    X = bsxfun(@times,Xp,BB.^(bta(ii,:)')');
    Y = bsxfun(@times,Yp,BB.^(bta(ii,:)')');
    
    score(ii) = KNN_obj(Rq, X, Y, obj_option, obj_para);
end




end

