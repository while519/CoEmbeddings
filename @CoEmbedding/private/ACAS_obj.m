function [score, Rqz, Rq] = ACAS_obj(R, Rz, obj_option, obj_para)
  
  switch lower(obj_option)
               
    case {'quant_integer', 'quant_int', 'quant_adapt',...
          'fast_quant_integer', 'fast_quant_int', 'fast_quant_adapt'}
      Rqz     =  ACAS_R(Rz, obj_option, obj_para);
      score   =  (sum(sum(( R-Rqz).^2))/size(Rqz,1)/size(Rqz,2))^0.5;
    case 'frobenius'
      score   = (sum(sum(( R-Rz).^2))/size(Rz,1)/size(Rz,2))^0.5;
      Rq  =R;
      Rqz =Rz;
    case 'alignment'
      score   = -trace(R*Rz')/sqrt(trace(R*R')*trace(Rz*Rz'));
      Rq  =R;
      Rqz =Rz;
    case {'log', 'log_likelihood', 'likelihood'}
      Rz= Rz +1e-15;
      score   = R.*log( Rz/(sum(Rz(:))) )/ sum(R(:));
      score   = -sum(sum(score));
      Rq  =R;
      Rqz =Rz;
    case {'kl', 'kl_divergence'}
      Rz= Rz +1e-15;
      A1 = diag(1./sum(R,2))*R;
      B1 = diag(1./sum(Rz,2))*Rz;
      score  = -sum(sum(A1.*log(B1)));  
      Rq  =R;
      Rqz =Rz;    
    otherwise
      error( [ 'unknown optimization objective function: ' obj_option ] );
  end


end