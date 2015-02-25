%%% CoEmbedding algorithm running script
clear classes
clear
clc
close all

path='/Users/yuwu/Documents/MATLAB/my_matlab/Relational_embedding/LED_release/datasets';
addpath(path);

%match={ 'ring','round','co-dot','SD1','SD2','SD3','SD4','SSD1', 'SSD2','SSD3', 'SSD4',...
%       'xor','2dnormals','Aggregation','Compound','circles', 'R15','rface','Pathbased',...
%       'CoTarget','CoIris','LetterRecognize','Newletter','ecoli','Codigit','digit04'}

load('SD3');
R=CalRelationXY(X,Y,'gaussian','mean');

figure
Plotcluster(X',lx,Y',ly);
axis equal


% %% CODE
% 
% opt  = {'optimization', 'ga', 'obj_option', 'quant_int', 'obj_para', {'bin', 20}};
% 
% I    = {'R', R', 'method', 'CODE', 'dim',2, 'model_optimisation', opt };
% o    = CoEmbedding(I);
% o    = training(o);
% 
% figure
% Plotcluster(o.X',lx,o.Y',ly);
% title('CODE Algorithm');
% 
% 
% %% ************************ACAS data plot********************************
% figure
% 
% opt  = {'optimization', 'grid', 'obj_option', 'quant_int', 'obj_para', {'bin', 10}};
% para_range = {'p', -1:3,  'alpha', 0:0.5:4, 'beta', 0:0.5:2};
% % I    = {'R', R, 'method', 'ACAS', 'dim',2, 'model_optimisation', opt,'parameter_range',para_range };
% I    = {'R', R, 'method', 'CA', 'dim',2};
% 
% o2    = CoEmbedding(I);
% o2    = training(o2);
% 
% figure
% Plotcluster(o2.Y',ly,o2.X',lx);
%              title('ACAS Algorithm');
             

% %% LEDF Model Optimization
% 
% opt  = {'optimisation', 'ga', 'obj_option', 'logic', 'obj_para', {'knn', [nan,nan]}};   
% 
% para_range = {'eta1',[-1,30],'eta2', [-1,30],'alpha',[0,3], 'beta', [-0.4,4] };
% 
% I    = {'R', R, 'method', 'LEDsvd', 'dim',2, 'model_optimisation', opt,'parameter_range',para_range };
% 
% 
% o    = CoEmbedding(I);
% o    = train(o);
% 
% h=figure;
% Plotcluster(o.X',lx,o.Y',ly);
% axis equal
% 
% return



opt  = {'optimisation', 'grid', 'obj_option', 'quant_knn', 'obj_para', {'knn', 10}};

para_range = {'eta1',[-1,9],'eta2', [-1,9],'alpha',[0.2,3], 'beta',[-0.5,1]};

I    = {'R', R, 'method', 'LEDeig', 'dim',2, 'model_optimisation', opt,'parameter_range',para_range };


o1    = CoEmbedding(I);
o1    = training(o1);

h=figure;
Plotcluster(o1.X',lx,o1.Y',ly);
axis equal


% %% CoEmbedding model evaluation
% 
% option  = {'evaluation', '3D', 'obj_option', 'quant_knn', 'obj_para', {'knn', 20}};
% 
% para_range = {'eta1', 1  , 'eta2', 2 , 'alpha', 0.2, 'beta', [0, 0.2], 'interval', 0.002 };
% 
% I    = {'R', R, 'method', 'LEDsvd', 'dim',2, 'model_optimisation', option,'parameter_range',para_range };
% 
% 
% obj_fun_eval = CoEmbedding(I);
% obj_fun_eval = Eval(obj_fun_eval);




