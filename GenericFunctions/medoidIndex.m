function index = medoidIndex(sequences, maxStretch) 
	

index = -1;
    lowestInertia = Inf;
	for i=1:length(sequences)
        tmpInertia = sumOfSquares(sequences{i},sequences, maxStretch);
        if (tmpInertia < lowestInertia)
            index = i;
            lowestInertia = tmpInertia;
        end
    end
    
    
    
    
    function sos = sumOfSquares(s,sequences, maxStretch)
    sos = 0.0;
    for j=1:length(sequences)
        dist = dtw(s,sequences{j}, maxStretch, 'squared');
        sos = sos + dist;
    end
end
    
    
    
end


