function [ result ] = rms( vector )

result = sqrt(mean(vector.^2));

end

