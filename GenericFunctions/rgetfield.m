function output = rgetfield(S, field)
    % Split the fieldname on "."
    parts = regexp(field, '\.', 'once', 'split');

    output = [S.(parts{1})];

    if numel(parts) > 1
        % If there are more parts, recursively get them
        output = rgetfield(output, parts{2});
    end
end