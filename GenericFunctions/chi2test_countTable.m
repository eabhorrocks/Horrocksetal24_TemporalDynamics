function [X2stat, p] = chi2test_countTable(countTable)

row_totals = sum(countTable,2);
col_totals = sum(countTable,1);

% calc expected vals
for irow = 1:size(countTable,1)
    for icol = 1:size(countTable,2)
        exp_table(irow,icol) = (row_totals(irow)*col_totals(icol))/sum(row_totals);
    end
end

tableSize = size(countTable);
df = (tableSize(1)-1)*(tableSize(2)-1);

observed = countTable(:); expected = exp_table(:);
X2stat = sum((observed-expected).^2./expected);
p = 1 - chi2cdf(X2stat, df);
