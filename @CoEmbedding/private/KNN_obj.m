function [score, Rqz, Rq] = KNN_obj(Rq, X, Y, obj_option, obj_para)

for ii = 1 : 2 : length(obj_para)-1
    switch lower(obj_para{ii})
        case 'bin'
            Rz = CalRelationXY(X, Y, 'gaussian', 'mean');
        case 'knn'
            if isscalar(obj_para{ii+1})
                knn_row = obj_para{ii+1};
                knn_col = obj_para{ii+1};
            else
                knn_row = obj_para{ii+1}(1);
                knn_col = obj_para{ii+1}(2);
            end
        otherwise
            error( [ 'unknown input argument: ' I{ii} ] );
    end
end

switch lower(obj_option)
    
    case 'logic'
        [M,N] = size(Rq);
        
        if ~isnan(knn_row)
            R1=zeros(M,N);
            % [~,idx]=sort(Dz,1);
            % per=bsxfun(@plus,idx(1:knn,:),M*(0:N-1));
            
            % [~,idx]=mink(Dz,knn);
            % per=bsxfun(@plus,idx,M*(0:N-1));
            [~,idx]=pdist2(X,Y,'euclidean','Smallest',knn_row);             % column
            per=bsxfun(@plus,idx,M*(0:N-1));
            R1(per(:))=1;
            R1 = sparse(R1);
        end
        
        if ~isnan(knn_col)
            % [~,idx]=sort(Dz',1);
            % per=bsxfun(@plus,idx(1:knn,:),N*(0:M-1));
            
            % [~,idx]=mink(Dz',knn);
            % per=bsxfun(@plus,idx,N*(0:M-1));
            
            R2=zeros(N,M);
            [~,idx]=pdist2(Y,X,'euclidean','Smallest',knn_col);
            per=bsxfun(@plus,idx(1:knn_col,:),N*(0:M-1));
            R2(per(:))=1;
            R2=sparse(R2)';
        end
        
        if isnan(knn_row) && ~isnan(knn_col)
            
            Rqz = R2.*Rq;
            
        elseif isnan(knn_col) && ~isnan(knn_row)
            
            Rqz = R1.*Rq;
            
        elseif ~isnan(knn_row) && ~isnan(knn_col)
            
            Rqz = R1.*R2.*Rq;
            
        else
            Rz = zeros(N,M);
            row_knn = sum(Rq,2);                                 % row vector of how much knn is needed in each column
            K = max(row_knn);
            [~,idx]=pdist2(Y,X,'euclidean','Smallest',K);
            for ii = 1:M
                Rz(idx(1:row_knn(ii),ii)+N*(ii-1)) = 1;
            end
            Rqz = Rz' ;
            
        end
        
        score=sum(abs(Rq(:)-Rqz(:)));
        
        
    case 'quant_knn'
        [M,N] = size(Rq);
        
        if ~isnan(knn_row)
            R1=zeros(M,N);
            % [~,idx]=sort(Dz,1);
            % per=bsxfun(@plus,idx(1:knn,:),M*(0:N-1));
            
            % [~,idx]=mink(Dz,knn);
            % per=bsxfun(@plus,idx,M*(0:N-1));
            [~,idx]=pdist2(X,Y,'euclidean','Smallest',knn_row);             % column
            per=bsxfun(@plus,idx,M*(0:N-1));
            R1(per(:))=1;
            R1 = sparse(R1);
        end
        
        if ~isnan(knn_col)
            % [~,idx]=sort(Dz',1);
            % per=bsxfun(@plus,idx(1:knn,:),N*(0:M-1));
            
            % [~,idx]=mink(Dz',knn);
            % per=bsxfun(@plus,idx,N*(0:M-1));
            
            R2=zeros(N,M);
            [~,idx]=pdist2(Y,X,'euclidean','Smallest',knn_col);
            per=bsxfun(@plus,idx(1:knn_col,:),N*(0:M-1));
            R2(per(:))=1;
            R2=sparse(R2)';
        end
        
        if isnan(knn_row)
            Rqz = R2.*Rq;
            
        elseif isnan(knn_col)
            Rqz = R1.*Rq;
        else
            Rqz = R1.*R2.*Rq;
        end
        
        score=sum(abs(Rq(:)-Rqz(:)));
        
    case 'quant_knn_abs'
        [M,N] = size(Rq);
        
        if ~isnan(knn_row)
            R1=zeros(M,N);
            % [~,idx]=sort(Dz,1);
            % per=bsxfun(@plus,idx(1:knn,:),M*(0:N-1));
            
            % [~,idx]=mink(Dz,knn);
            % per=bsxfun(@plus,idx,M*(0:N-1));
            [~,idx]=pdist2(X,Y,'euclidean','Smallest',knn_row);             % column
            per=bsxfun(@plus,idx,M*(0:N-1));
            R1(per(:))=1;
            R1 = sparse(R1);
        end
        
        if ~isnan(knn_col)
            % [~,idx]=sort(Dz',1);
            % per=bsxfun(@plus,idx(1:knn,:),N*(0:M-1));
            
            % [~,idx]=mink(Dz',knn);
            % per=bsxfun(@plus,idx,N*(0:M-1));
            
            R2=zeros(N,M);
            [~,idx]=pdist2(Y,X,'euclidean','Smallest',knn_col);
            per=bsxfun(@plus,idx(1:knn_col,:),N*(0:M-1));
            R2(per(:))=1;
            R2=sparse(R2)';
        end
        
        if isnan(knn_row)
            Rqz = R2;
            
        elseif isnan(knn_col)
            Rqz = R1;
        else
            Rqz = R1.*R2;
        end
        
        score=sum(abs(Rq(:)-Rqz(:)));
        
    case {'quant_integer', 'quant_int', 'quant_adapt',...
            'fast_quant_integer', 'fast_quant_int', 'fast_quant_adapt'}
        Rqz     =  ACAS_R(Rz, obj_option, obj_para);
        score   =  (sum(sum(( Rq-Rqz).^2))/size(Rqz,1)/size(Rqz,2))^0.5;
    case 'frobenius'
        score   = (sum(sum(( Rq-Rz).^2))/size(Rz,1)/size(Rz,2))^0.5;
        Rqz =Rz;
    case 'alignment'
        score   = -trace(Rq*Rz')/sqrt(trace(Rq*Rq')*trace(Rz*Rz'));
        Rqz =Rz;
    case {'log', 'log_likelihood', 'likelihood'}
        Rz= Rz +1e-15;
        score   = Rq.*log( Rz/(sum(Rz(:))) )/ sum(Rq(:));
        score   = -sum(sum(score));
        Rqz =Rz;
    case {'kl', 'kl_divergence'}
        Rz= Rz +1e-15;
        A1 = diag(1./sum(Rq,2))*Rq;
        B1 = diag(1./sum(Rz,2))*Rz;
        score  = -sum(sum(A1.*log(B1)));
        Rqz =Rz;
    otherwise
        error( [ 'unknown optimization objective function: ' obj_option ] );
end


end