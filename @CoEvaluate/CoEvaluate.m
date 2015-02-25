classdef CoEvaluate
    
    properties
        ScoreType = 'Within_Group';
        X
        Y
        lx
        ly
        knn = 7;
        k_fold = 10;
        BarMat
    end % properties
    
    
    properties (Dependent = true, SetAccess = private )
        SW                       % Within group classification score
        SB                       % Between group classification score
        SS                       % Cross validation accuracy rate
        SM                       % score in [W, B, S] row vector format
    end % properties
    methods
        function ce = CoEvaluate(X, lx, Y, ly)
            if nargin > 0   % Support calling with 0 argument
                ce.X = X;
                ce.Y = Y;
                ce.lx = lx;
                ce.ly = ly;
            end
        end % CoEvaluate
        
        function disp(obj)
            fprintf(1, 'Within Group score:  %6.4f\nBetween Group score:  %6.4f\nCross Validate score:  %6.4f\n',...
                obj.SW, obj.SB, obj.SS);
        end % disp W,S,B score
        
        function bar(obj,varargin)
            figure
            ax = axes;
            bar(obj.BarMat,varargin{:});
            legend('Origin','LEDF', 'ACAS', 'CODE');
            ax.XTickLabel = {'W-score','B-score','S-score'};
            ax.YLim = [0,1];
        end % bar chart
        
        function obj = set.ScoreType(obj,score)
            score = lower(score);
            if ~(strcmp(score, 'within_group') || ...
                    strcmp(score, 'between_group') || ...
                    strcmp(score, 'cross_validate'))
                error('score type must be W, S, W');
            end
            obj.ScoreType = score;
        end % set.ScoreType
        
        function sw = get.SW(obj)
            % This score evaluate the class separability for X and Y group that it
            % measures the how much the same class label data in X
            % and Y group stay close
            
            % for s in X, measure the within group score
            
            Wx = 0;
            
            for ii = 1 : length(obj.lx)
                s = obj.X(ii,:);
                ls = obj.lx(ii);
                
                fsetx = find(obj.lx == ls);      % the friendly sets for s including s itself
                
                t = 0;
                
                for kk = 1: obj.knn
                    [~,IDX] = pdist2(obj.X , s , 'euclidean', 'Smallest', kk+1);
                    Lia = ismember(IDX,fsetx);
                    t = t+(sum(Lia)-1)/kk;
                end
                
                Wx = Wx+t/obj.knn;
                
            end
            
            % for s in Y, measure the within group score
            Wy = 0;
            for jj = 1 : length(obj.ly)
                s = obj.Y(jj,:);
                ls = obj.ly(jj);
                
                fsety = find(obj.ly == ls);     % friendly sets to s in Y group
                
                t=0;
                
                
                for kk = 1: obj.knn
                    [~,IDY] = pdist2(obj.Y , s, 'euclidean','Smallest', kk+1);
                    Lia = ismember(IDY, fsety);
                    t = t+(sum(Lia)-1)/kk;
                end
                
                Wy = Wy+t/obj.knn;
            end
            
            str=['Wx =' num2str(Wx/length(obj.lx)) ' Wy =' num2str(Wy/length(obj.ly))];
            disp(str)
            sw =(Wx+Wy)/(length(obj.lx)+length(obj.ly));
        end % get.sw
        
        function sb = get.SB(obj)
            % This score evaluate the inter-group separability for co-class in
            % X and Y group that it measures the how much the inter-group co-class pair
            % stay close one certain co-class
            
            % Identify the co-clustering data label in X group
            Ix = find(ismember(obj.lx, obj.ly) == 1);
            
            if isempty(Ix)
                warning('no co-clustering data detect in this sample');
                sb = 0;
                return
            end
            
            B1x = 0;
            
            % for s in one co-cluster in X group
            for ii = 1 : length(Ix)
                s = obj.X(Ix(ii),:);
                ls = obj.lx(Ix(ii));
                
                fsets = find(obj.ly == ls);     % friend sets to s in Y group
                
                t=0;
                
                for kk = 1 : obj.knn
                    [~, IDX] = pdist2(obj.Y , s , 'euclidean', 'Smallest', kk);
                    Lia = ismember(IDX , fsets);
                    t = t+(sum(Lia))/kk;
                end
                
                B1x = B1x+t/obj.knn;
                
            end
            
            
            Iy = find(ismember(obj.ly, obj.lx) == 1);
            B1y = 0;
            
            % for s in one co-cluster in Y group
            for jj = 1 : length(Iy)
                s = obj.Y(Iy(jj),:);
                ls = obj.ly(Iy(jj));
                
                fsets = find(obj.lx == ls);     % friend sets to s in X group
                
                t = 0;
                
                for kk = 1 : obj.knn
                    [~, IDY] = pdist2(obj.X , s , 'euclidean', 'Smallest', kk);
                    Lia = ismember(IDY ,fsets);
                    t = t+(sum(Lia))/kk;
                end
                
                B1y = B1y+t/obj.knn;
                
            end
            
            str = ['B1x = ' num2str(B1x/length(Ix)) ' B1y = ' num2str(B1y/length(Iy))];
            disp(str)
            
            sb = (B1x+B1y)/(length(Iy)+length(Ix));
            
        end % get.sb
        
        function ss = get.SS(obj)
            rng('default');
            % cross validate using knn method
            mdl = fitcknn([obj.X; obj.Y], [obj.lx; obj.ly], 'NumNeighbors', obj.knn);
            cvmdl = crossval(mdl,'KFold', obj.k_fold);
            ss = 1-kfoldLoss(cvmdl);
        end % get.ss
        
        function sm = get.SM(obj)
            sm = [obj.SW, obj.SB, obj.SS];
        end
        
    end % methods
end % classdef