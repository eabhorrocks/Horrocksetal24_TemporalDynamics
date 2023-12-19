function cosTheta = calcCosineTheta(v1,v2)

cosTheta = dot(v1,v2)./(norm(v1)*norm(v2));