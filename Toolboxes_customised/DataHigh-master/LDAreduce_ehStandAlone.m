function [newD, C] = LDAreduce_ehStandAlone(D,dims)
%LDAREDUCE Internal function for LDA
%   LDAREDUCE(D,DIMS) returns a structure of the same form as D, except the
%   data has been reduced with LDA. All conditions and trials are
%   considered together to get the best joint reduction.
    conds = length(unique({D.condition}));
%     if dims > conds
%         return;
%     end
    
   % try
        [newD lda_eigs] = lda_engineDH(D,dims);
        C = lda_eigs(:,1:dims);
        params.d = mean([D.data],2);
    %catch err
    %    fprintf(['\n\nLDA failed.  Check to make sure you have more than\n' ...
%             'condition and more than one trial per condition.\n\n']);
%    end

end