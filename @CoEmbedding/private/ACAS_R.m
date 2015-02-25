function Rq= ACAS_R(R, type, para)
bin = 10;

for ii = 1 : 2 : length(para)-1
    switch lower(para{ii})
        case 'bin'
            bin	= para{ii+1};
        case 'knn'
            if isscalar(para{ii+1})
                knn_row = para{ii+1};
                knn_col = para{ii+1};
            else
                knn_row = para{ii+1}(1);
                knn_col = para{ii+1}(2);
            end
        otherwise
            error( [ 'unknown input argument: ' I{ii} ] );
    end
end

switch lower(type)
    
    case 'logic'
        Rq = sparse(R);
        
    case {'quant_knn','quant_knn_abs'}
        
        [M,N]=size(R);
        if ~isnan(knn_row)
            
            R1=zeros(M,N);
            
            [~,idx]=sort(R,1,'descend');                    % column
            per=bsxfun(@plus,idx(1:knn_row,:),M*(0:N-1));
            R1(per(:))=1;
            R1 = sparse(R1);
        end
        
        if ~isnan(knn_col)
            R2=zeros(N,M);
            [~,idx]=sort(R',1,'descend');                   % row
            per=bsxfun(@plus,idx(1:knn_col,:),N*(0:M-1));
            R2(per(:))=1;
            R2=sparse(R2)';
        end
        
        if isnan(knn_row)
            Rq = R2;
        elseif isnan(knn_col)
            Rq = R1;
        else
            Rq=R1.*R2;
        end
        
    case {'quant_integer','quant_int'}
        bound = quantile(R(:),bin);
        Rq=zeros(size(R));
        for ii=1:length(bound)
            if ii==1
                ind = find(R<=bound(1));
            elseif ii==length(bound)
                ind = find(R>bound(end));
            else
                ind = find( R>bound(ii-1) & R<=bound(ii) );
            end
            Rq(ind)=ii;
        end
    case {'fast_quant_integer','fast_quant_int'}
        tt = R(:);
        [tt,N] = sort(tt);
        if length(tt)>5000
            Nstep = round(length(tt)/5000)+1;
        else
            Nstep=1;
        end
        bound = quantile(tt(1:Nstep:length(tt)),bin);
        Rq=zeros(size(R));
        for ii=1:length(bound)
            if ii==1
                ind = find(R<=bound(1));
            elseif ii==length(bound)
                ind = find(R>bound(end));
            else
                ind = find( R>bound(ii-1) & R<=bound(ii) );
            end
            Rq(ind)=ii;
        end
    case 'quant_adapt'
        R = R/(sum(R(:)));
        bound = quantile(R(:),bin);
        Rq=zeros(size(R));
        for ii=1:length(bound)
            if ii==1
                ind = find(R<=bound(1));
                Rq(ind) = (min(R(:)) + bound(ii))/2;
            elseif ii==length(bound)
                ind = find(R>bound(end));
                Rq(ind) = (max(R(:)) + bound(ii))/2;
            else
                ind = find( R>bound(ii-1) & R<=bound(ii) );
                Rq(ind) = (bound(ii-1) + bound(ii))/2;
            end
        end
    case 'fast_quant_adapt'
        R = R/(sum(R(:)));
        tt = R(:);
        [tt,N] = sort(tt);
        if length(tt)>5000
            Nstep = round(length(tt)/5000)+1;
        else
            Nstep=1;
        end
        bound = quantile(tt(1:Nstep:length(tt)),bin);
        Rq=zeros(size(R));
        for ii=1:length(bound)
            if ii==1
                ind = find(R<=bound(1));
                Rq(ind) = (min(R(:)) + bound(ii))/2;
            elseif ii==length(bound)
                ind = find(R>bound(end));
                Rq(ind) = (max(R(:)) + bound(ii))/2;
            else
                ind = find( R>bound(ii-1) & R<=bound(ii) );
                Rq(ind) = (bound(ii-1) + bound(ii))/2;
            end
        end
    otherwise
        error( [ 'unknown method for processing ACAS optimisation objective: ',  type{ii} ] );
end

end

